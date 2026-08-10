#include "config.h"

#include <algorithm>
#include <charconv>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <stdexcept>
#include <tuple>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <string>
#include <string_view>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <curl/curl.h>

#include "date.h"
#include "ephemeris.h"
#include "exceptions.h"
#include "files.h"
#include "horizons.h"
#include "logging.h"
#include "orbital_elements.h"
#include "util/string.h"

using std::endl;
using std::string;
using std::string_view;

namespace fs = boost::filesystem;

namespace sbsearch
{
    Horizons::Horizons(const MovingTarget target,
                       string_view center,
                       const Date start_date,
                       const Date stop_date,
                       string_view step_size,
                       const bool cache)
        : target_(target),
          center_(center),
          start_date_(start_date),
          stop_date_(stop_date),
          step_size_(step_size),
          cache_(cache) {}

    string_view Horizons::command()
    {
        format_command();
        return command_;
    }

    string_view Horizons::parameters()
    {
        format_query();
        return parameters_;
    }

    string Horizons::format_command(const MovingTarget &target, const double mjd)
    {
        // for major body search, we're done!
        if (!target.small_body())
            return Horizons::major_body_command(target);

        if (target.orbit())
            return Horizons::orbit_command(target);

        return Horizons::small_body_command(target, mjd);
    }

    string Horizons::major_body_command(const MovingTarget &target)
    {
        return Horizons::command(target.designation());
    }

    string Horizons::small_body_command(const MovingTarget &target, const double mjd)
    {
        const string designation(target.designation());

        // Temporary comet designation?  C/2001 Q4, P/2003 CC22
        bool temporary_comet = (designation.find_first_of("CPDI") == 0) && (designation[1] == '/');

        // Permanent comet or interstellar object designation? 2P, 2I, 100P.
        // Check that there is a parsable number without whitespace before the
        // letter.
        bool comet = temporary_comet;
        if (!temporary_comet)
        {
            size_t letter = designation.find_first_of("DPI");
            if ((letter != string::npos) && !std::isspace(designation[letter - 1]))
                comet = std::all_of(designation.begin(), std::next(designation.begin(), letter), [](char c)
                                    { return std::isdigit(c); });
        }

        return comet ? Horizons::comet_command(target, mjd) : Horizons::asteroid_command(target);
    }

    string Horizons::asteroid_command(const MovingTarget &target)
    {
        string designation(target.designation());
        string s("");
        bool no_digits = designation.find_first_of("0123456789") == string::npos;
        bool all_digits = designation.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ") == string::npos;
        if (no_digits || all_digits)
            // probably an asteroid name or number
            s = designation;
        else
            // otherwise, assume it is a temporary ID
            s = "DES=" + designation;

        return Horizons::command(s + ";");
    }

    string Horizons::comet_command(const MovingTarget &target, const double mjd)
    {
        std::stringstream s;

        if (mjd > 0)
        {
            int jd = (int)(mjd + 2400000.5);
            s << "DES=" << target.designation() << ";NOFRAG;CAP<" << std::to_string(jd) << ";";
        }
        else
            s << "DES=" << target.designation() << ";NOFRAG;CAP;";

        return Horizons::command(s.str());
    }

    string Horizons::orbit_command(const MovingTarget &target)
    {
        /*
        OBJECT	 	Name of user input object
        EPOCH	 	Julian Day number (JDTDB) of osculating elements
        ECLIP	 	Reference ecliptic frame of elements: J2000 or B1950. J2000 assumes the IAU76/80 J2000 obliquity of 84381.448 arcsec relative to the ICRF reference frame. B1950 assumes FK4/B1950 obliquity of 84404.8362512 arcsec.
        EC	 	Eccentricity
        QR	au	Perihelion distance (see note above)
        TP	 	Perihelion Julian Day number (see note above)
        OM	deg	Longitude of ascending node wrt ecliptic
        W	deg	Argument of perihelion wrt ecliptic
        IN	deg	Inclination wrt ecliptic
        MA	deg	Mean anomaly (see note above)
        A	au	Semi-major axis (see note above)
        N	deg/d	Mean motion (see note above)
        */
        auto orbit = target.orbit().value();
        std::stringstream s;
        constexpr auto precision{std::numeric_limits<long double>::digits10 + 1};
        s << std::setprecision(precision)
          << Horizons::command(";") << "\n"
          << "OBJECT='" << target.designation() << "'\n"
          << "EPOCH=" << orbit.epoch << "\n"
          << "ECLIP=J2000" << "\n"
          << "EC=" << orbit.ec << "\n"
          << "OM=" << orbit.om << "\n"
          << "W=" << orbit.w << "\n"
          << "IN=" << orbit.in;

        // Could apply some logic here to make sure a consistent set of elements
        // are provided: {TP, QR}, {MA, A} or {MA,N} may be specified.  Horizons
        // manual says, "Note that if you specify elements with MA, {TP, QR}
        // will be computed from them. The program always uses TP and QR
        // internally."
        if (orbit.qr)
            s << "\n"
              << "QR=" << orbit.qr.value();

        if (orbit.Tp)
            s << "\n"
              << "TP=" << orbit.Tp.value();

        if (orbit.ma)
            s << "\n"
              << "MA=" << orbit.ma.value();

        if (orbit.a)
            s << "\n"
              << "A=" << orbit.a.value();

        if (orbit.n)
            s << "\n"
              << "N=" << orbit.n.value();

        return s.str();
    }

    void Horizons::format_command()
    {
        command_ = format_command(target_, stop_date_.mjd());
    }

    string Horizons::format_query(string_view command,
                                  string_view center,
                                  const Date &start_date,
                                  const Date &stop_date,
                                  string_view step_size)
    {
        std::stringstream stream;
        stream << "!$$SOF\n"
               << "MAKE_EPHEM=YES\n"
               << command << "\n"
               << "EPHEM_TYPE=OBSERVER\n"
               << "CENTER='" << center << "'\n"
               << "START_TIME='" << start_date.iso() << "'\n"
               << "STOP_TIME='" << stop_date.iso() << "'\n"
               << "STEP_SIZE='" << step_size << "'\n"
               << "QUANTITIES='1,9,19,20,23,24,27,37,41,47'\n"
               << "REF_SYSTEM='ICRF'\n"
               << "CAL_FORMAT='JD'\n"
               << "CAL_TYPE='M'\n"
               << "TIME_DIGITS='MINUTES'\n"
               << "ANG_FORMAT='DEG'\n"
               << "APPARENT='AIRLESS'\n"
               << "RANGE_UNITS='AU'\n"
               << "SUPPRESS_RANGE_RATE='NO'\n"
               << "SKIP_DAYLT='NO'\n"
               << "SOLAR_ELONG='0,180'\n"
               << "EXTRA_PREC='YES'\n"
               << "R_T_S_ONLY='NO'\n"
               << "CSV_FORMAT='YES'\n"
               << "OBJ_DATA='YES'\n";

        return stream.str();
    }

    void Horizons::format_query()
    {
        // string_view cmd = command();
        parameters_ = format_query(command(), center_, start_date_, stop_date_, step_size_);
    }

    string Horizons::query(string_view parameters, const bool cache)
    {
        const fs::path fn = generate_cache_file_name(parameters);
        if (cache && fs::exists(fn))
        {
            Logger::info() << "Reading Horizons cache: " << fn << endl;
            return read_file(fn.string());
        }

        Logger::debug() << "Query Horizons with parameters\n"
                        << parameters << endl;

        string table;
        try
        {
            char error_message[CURL_ERROR_SIZE];

            CURL *handle = curl_easy_init();
            curl_easy_setopt(handle, CURLOPT_URL, "https://ssd.jpl.nasa.gov/api/horizons_file.api");
            curl_easy_setopt(handle, CURLOPT_USERAGENT, "SBSearch/" SBSEARCH_VERSION);
            curl_easy_setopt(handle, CURLOPT_WRITEFUNCTION, write_http_string_data);
            curl_easy_setopt(handle, CURLOPT_WRITEDATA, (void *)&table);
            curl_easy_setopt(handle, CURLOPT_ERRORBUFFER, error_message);

            curl_mime *multipart = curl_mime_init(handle);
            curl_mimepart *part = curl_mime_addpart(multipart);
            curl_mime_name(part, "format");
            curl_mime_data(part, "text", CURL_ZERO_TERMINATED);
            part = curl_mime_addpart(multipart);
            curl_mime_name(part, "input");
            curl_mime_data(part, parameters.data(), parameters.size());
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

                throw HorizonsError(user_message);
            }

            curl_easy_cleanup(handle);

            const string api_version = "API VERSION: 1.0";
            if (table.find(api_version) == string::npos)
                throw HorizonsError("Expected version string: " + api_version);
            // test for an ephemeris before writing to the cache
            int pos = table.find("$$SOE\n");
            if (pos == string::npos)
                throw HorizonsError("Start of ephemeris string ($$SOE) not found in data table.");

            pos = table.find("$$EOE\n");
            if (pos == string::npos)
                throw HorizonsError("End of ephemeris string ($$EOE) not found in data table.");
        }
        catch (std::exception &e)
        {
            std::cerr << parameters << "\n\n"
                      << table;
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

    Ephemeris::Data Horizons::parse(string_view table)
    {
        Ephemeris::Data data;

        // find the start of the data
        int data_start = table.find("$$SOE\n");
        if (data_start == string::npos)
            throw HorizonsError("Start of ephemeris string ($$SOE) not found in data table.");
        data_start += 6;

        const int data_end = table.find("$$EOE\n");
        if (data_end == string::npos)
            throw HorizonsError("End of ephemeris string ($$EOE) not found in data table.");

        // Find the period and time of perihelion for T-Tp calculations
        double period = 0, Tp = 0;
        int i;
        if ((i = table.find("TP=")) != string::npos)
        {
            auto [ptr, ec] = std::from_chars(table.data() + i + 3,
                                             table.data() + table.size() - i - 3,
                                             Tp);
            if (ec != std::errc())
                throw std::invalid_argument("Cannot parse Tp value as double");
        }

        if ((i = table.find("PER=")) != string::npos)
        {
            // period = std::stod(table.substr(i + 4)) * 365.25;
            auto [ptr, ec] = std::from_chars(table.data() + i + 4,
                                             table.data() + table.size() - i - 4,
                                             period);
            if (ec != std::errc())
                throw std::invalid_argument("Cannot parse PER value as double");

            period *= 365.25;
        }

        // get column names, 3 lines before the first data line
        int column_names_start = data_start;
        for (int i = 0; i < 3; i++)
            column_names_start = table.rfind("\n", column_names_start - 2) + 1;

        vector<string> column_names = util::split(
            table.substr(
                column_names_start,
                table.find("\n", column_names_start) - column_names_start),
            ',');

        // Remove whitespace and underscores from column names.  The latter can vary
        // if the column length varies (e.g., when the extended precision flag is
        // enabled)

        auto whitespace_or_underscore = [](int c)
        { return std::isspace(c) || (c == '_'); };

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
            "r", "delta", "S-T-O", "S-O-T", "TruAnom", "PsAng", "PsAMV", "Skymotion",
            "SkymotPA"};
        for (string name : data_names)
        {
            auto i = std::find(column_names.begin(), column_names.end(), name);
            if (i == column_names.end())
                throw HorizonsError("Column " + name + " not found in data table.");
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
            throw HorizonsError("No magnitude columns found, searched for T-mag, N-mag, and APmag.");

        // iterate over rows of data
        int row_start = data_start;

        // Convert cell to double, but if n.a., return 0
        auto celltod = [](string_view s)
        { return (s.find("n.a.") == string::npos) ? util::svtod(s) : 0; };

        while (true)
        {
            int line_length = table.find("\n", row_start) - row_start;
            string_view line = table.substr(row_start, line_length);

            vector<string> row = util::split(line, ',');

            double vmag = 99;
            for (int i : magnitude_column_indices)
            {
                string value = row[i];
                if (value.find("n.a.") != string::npos)
                    continue;
                vmag = (std::stod(value) < vmag) ? std::stod(value) : vmag;
            }

            double jd = std::stod(row[columns["DateJDUT"]]);

            // Horizons's Theta is clockwise from east.  Change to
            // counter-clockwise from north.
            data.push_back({jd - 2400000.5,
                            std::fmod(jd - Tp, period),
                            std::stod(row[columns["R.A.(ICRF)"]]),
                            std::stod(row[columns["DEC(ICRF)"]]),
                            std::stod(row[columns["Skymotion"]]),
                            std::stod(row[columns["SkymotPA"]]),
                            celltod(row[columns["SMAA3sig"]]),
                            celltod(row[columns["SMIA3sig"]]),
                            360 - celltod(row[columns["Theta"]]) - 90,
                            std::stod(row[columns["r"]]),
                            std::stod(row[columns["delta"]]),
                            std::stod(row[columns["S-T-O"]]),
                            std::stod(row[columns["S-O-T"]]),
                            celltod(row[columns["TruAnom"]]),
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
        Logger::info() << "Querying Horizons for ephemeris: " << target_.to_string()
                       << " from " << start_date_.iso()
                       << " to " << stop_date_.iso()
                       << " with step size " << step_size_
                       << endl;
        query();
        parse();
        return data_;
    }

    string Horizons::command(string_view s)
    {
        std::stringstream stream;
        stream << "COMMAND='" << s << "'";
        return stream.str();
    };
}
