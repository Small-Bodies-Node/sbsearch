#include "config.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <stdexcept>
#include <tuple>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <curl/curl.h>

#include "date.h"
#include "ephemeris.h"
#include "exceptions.h"
#include "files.h"
#include "horizons.h"
#include "logging.h"

using std::string;
namespace fs = boost::filesystem;

namespace sbsearch
{
    Horizons::Horizons(const MovingTarget target,
                       const string center,
                       const Date start_date,
                       const Date stop_date,
                       const string time_step,
                       const bool cache)
    {
        target_ = target;
        center_ = center;
        start_date_ = start_date;
        stop_date_ = stop_date;
        time_step_ = time_step;
        cache_ = cache;
    }

    string Horizons::format_command(const string designation,
                                    const bool small_body,
                                    const double mjd)
    {
        // for major body search, we're done!
        if (!small_body)
            return designation;

        // Temporary comet designation?  C/2001 Q4, P/2003 CC22
        bool temporary_comet = (designation.find_first_of("CPDI") == 0) & (designation[1] == '/');

        // Permanent comet or interstellar object designation? 2P, 2I, 100P.
        // Check that there is a parsable number without whitespace before the
        // letter.
        bool comet = temporary_comet;
        if (!temporary_comet)
        {
            size_t letter = designation.find_first_of("DPI");
            if ((letter != string::npos) & !std::isspace(designation[letter - 1]))
            {
                try
                {
                    std::stoul(designation.substr(0, letter));
                    comet = true;
                }
                catch (std::exception &e)
                {
                    comet = false;
                }
            }
        }

        if (comet)
        {
            if (mjd > 0)
            {
                int jd = (int)(mjd + 2400000.5);
                return "DES=" + designation + ";NOFRAG;CAP<" + std::to_string(jd) + ";";
            }
            else
                return "DES=" + designation + ";NOFRAG;CAP;";
        }

        if (designation.find_first_of("0123456789") == string::npos)
            // probably an asteroid name, just append a semicolon
            return designation + ";";

        if (designation.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ") == string::npos)
            // probably a permanent asteroid number, just append a semicolon
            return designation + ";";

        // otherwise, assume it is a temporary ID
        return "DES=" + designation + ";";
    }

    void Horizons::format_command()
    {
        command_ = format_command(target_.designation(), target_.small_body(), start_date_.mjd());
    }

    string Horizons::format_query(const string command,
                                  const string center,
                                  const Date start_date,
                                  const Date stop_date,
                                  const string time_step)
    {
        char parameters[1024];
        sprintf(parameters, R"(
!$$SOF
MAKE_EPHEM=YES
COMMAND='%s'
EPHEM_TYPE=OBSERVER
CENTER='%s'
START_TIME='%s'
STOP_TIME='%s'
STEP_SIZE='%s'
QUANTITIES='1,9,19,20,23,24,27,37,41'
REF_SYSTEM='ICRF'
CAL_FORMAT='JD'
CAL_TYPE='M'
TIME_DIGITS='MINUTES'
ANG_FORMAT='DEG'
APPARENT='AIRLESS'
RANGE_UNITS='AU'
SUPPRESS_RANGE_RATE='NO'
SKIP_DAYLT='NO'
SOLAR_ELONG='0,180'
EXTRA_PREC='YES'
R_T_S_ONLY='NO'
CSV_FORMAT='YES'
OBJ_DATA='YES'
)",
                command.c_str(),
                center.c_str(),
                start_date.iso().c_str(),
                stop_date.iso().c_str(),
                time_step.c_str());

        return string(parameters);
    }

    void Horizons::format_query()
    {
        parameters_ = format_query(command(), center_, start_date_, stop_date_, time_step_);
    }

    string Horizons::query(const string parameters, const bool cache)
    {
        const fs::path fn = generate_cache_file_name(parameters);
        if (cache & fs::exists(fn))
            return read_file(fn.string());

        string table;
        try
        {
            char error_message[CURL_ERROR_SIZE];

            CURL *handle = curl_easy_init();
            curl_easy_setopt(handle, CURLOPT_URL, "https://ssd.jpl.nasa.gov/api/horizons_file.api");
            curl_easy_setopt(handle, CURLOPT_WRITEFUNCTION, write_http_string_data);
            curl_easy_setopt(handle, CURLOPT_WRITEDATA, (void *)&table);
            curl_easy_setopt(handle, CURLOPT_ERRORBUFFER, error_message);

            curl_mime *multipart = curl_mime_init(handle);
            curl_mimepart *part = curl_mime_addpart(multipart);
            curl_mime_name(part, "format");
            curl_mime_data(part, "text", CURL_ZERO_TERMINATED);
            part = curl_mime_addpart(multipart);
            curl_mime_name(part, "input");
            curl_mime_data(part, parameters.c_str(), CURL_ZERO_TERMINATED);
            curl_easy_setopt(handle, CURLOPT_MIMEPOST, multipart);

            // Consider CURLOPT_VERBOSE and CURLOPT_DEBUGFUNCTION to better debug and trace why errors happen.
            // if (verbose)
            //     curl_easy_setopt(handle, CURLOPT_VERBOSE, 1);

            CURLcode code = curl_easy_perform(handle);

            if (code != CURLE_OK)
            {
                char user_message[CURL_ERROR_SIZE + 256];
                size_t len = std::strlen(error_message);
                if (len)
                    sprintf(user_message, "Error fetching JPL Horizons data: %s", error_message);
                else
                    sprintf(user_message, "%s", curl_easy_strerror(code));

                throw std::runtime_error(user_message);
            }

            curl_easy_cleanup(handle);

            const string api_version = "API VERSION: 1.0";
            if (table.find(api_version) == string::npos)
                throw std::runtime_error("Unexpected Horizons response version.  Expected: " + api_version);
        }
        catch (std::exception &e)
        {
            std::cerr << table;
            throw;
        }

        if (cache)
            write_to_cache(fn, table);

        return table;
    }

    void Horizons::query()
    {
        table_ = query(parameters(), cache_);
    }

    Ephemeris::Data Horizons::parse(const string &table)
    {
        Ephemeris::Data data;

        // find the start of the data
        int data_start = table.find("$$SOE\n");
        if (data_start == string::npos)
            throw std::runtime_error("Start of ephemeris string ($$SOE) not found in data table.");
        data_start += 6;

        const int data_end = table.find("$$EOE\n");
        if (data_end == string::npos)
            throw std::runtime_error("End of ephemeris string ($$EOE) not found in data table.");

        // Find the period and time of perihelion for T-Tp calculations
        double period = 0, Tp = 0;
        int i;
        if ((i = table.find("TP=")) != string::npos)
            Tp = std::stod(table.substr(i + 3));
        if ((i = table.find("PER=")) != string::npos)
            period = std::stod(table.substr(i + 4)) * 365.25;

        // get column names, 3 lines before the first data line
        int column_names_start = data_start;
        for (int i = 0; i < 3; i++)
            column_names_start = table.rfind("\n", column_names_start - 2) + 1;

        vector<string> column_names = sbsearch::split(
            table.substr(
                column_names_start,
                table.find("\n", column_names_start) - column_names_start),
            ',');

        // Remove whitespace and underscores from column names.  The latter can vary
        // if the column length varies (e.g., when the extended precision flag is
        // enabled)

        auto whitespace_or_underscore = [](int c)
        { return std::isspace(c) | (c == '_'); };

        auto remove_whitespace_and_underscores = [whitespace_or_underscore](string s)
        {
            auto end = std::remove_if(s.begin(), s.end(), whitespace_or_underscore);
            return s.substr(0, end - s.begin());
        };

        std::transform(column_names.begin(), column_names.end(), column_names.begin(), remove_whitespace_and_underscores);

        // Store the indices of the columns we are interested in.
        std::map<string, int> columns;
        const vector<string> data_names{
            "DateJDUT", "R.A.(ICRF)", "DEC(ICRF)", "SMAA3sig", "SMIA3sig", "Theta",
            "r", "delta", "S-T-O", "S-O-T", "TruAnom", "PsAng", "PsAMV"};
        for (string name : data_names)
        {
            auto i = std::find(column_names.begin(), column_names.end(), name);
            if (i == column_names.end())
                throw std::runtime_error("Column " + name + " not found in data table.");
            columns[name] = i - column_names.begin();
        }

        // find any magnitude columns
        vector<int> magnitude_column_indices;
        for (string name : {"T-mag", "N-mag", "APmag"})
        {
            auto i = std::find(column_names.begin(), column_names.end(), name);
            if (i == column_names.end())
                continue;
            magnitude_column_indices.push_back(i - column_names.begin());
        }
        if (magnitude_column_indices.size() == 0)
            throw std::runtime_error("No magnitude columns found, searched for T-mag, N-mag, and APmag.");

        // iterate over rows of data
        int row_start = data_start;

        // Convert cell to double, but if n.a., return 0
        auto celltod = [](const string &s)
        { return (s.find("n.a.") == string::npos) ? std::stod(s) : 0; };

        while (true)
        {
            int line_length = table.find("\n", row_start) - row_start;
            string line = table.substr(row_start, line_length);

            vector<string> row = sbsearch::split(line, ',');

            double vmag = 99;
            for (int i : magnitude_column_indices)
            {
                string value = row[i];
                if (value.find("n.a.") != string::npos)
                    continue;
                vmag = (std::stod(value) < vmag) ? std::stod(value) : vmag;
            }

            double jd = std::stod(row[columns["DateJDUT"]]);
            data.push_back({jd - 2400000.5,
                            std::fmod(jd - Tp, period),
                            std::stod(row[columns["R.A.(ICRF)"]]),
                            std::stod(row[columns["DEC(ICRF)"]]),
                            celltod(row[columns["SMAA3sig"]]),
                            celltod(row[columns["SMIA3sig"]]),
                            celltod(row[columns["Theta"]]),
                            std::stod(row[columns["r"]]),
                            std::stod(row[columns["delta"]]),
                            std::stod(row[columns["S-T-O"]]),
                            std::stod(row[columns["S-O-T"]]),
                            std::stod(row[columns["TruAnom"]]),
                            std::fmod(std::stod(row[columns["PsAng"]]) + 180., 360),
                            std::fmod(std::stod(row[columns["PsAMV"]]) + 180., 360),
                            vmag});

            row_start += line_length + 1;
            if (row_start >= data_end)
                break;
        }

        return data;
    }

    void Horizons::parse()
    {
        data_ = parse(table_);
    }

    Ephemeris::Data Horizons::get_ephemeris_data()
    {
        query();
        parse();
        return data_;
    }
}
