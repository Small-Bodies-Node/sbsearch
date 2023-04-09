#include <fstream>
#include <iostream>
#include <cstdio>
#include <streambuf>
#include <string>
#include <vector>
#include <boost/program_options.hpp>
#include <curl/curl.h>

#include "ephemeris.h"
#include "sbs-cli.h"
#include "sofa/sofa.h"

using std::cerr;
using std::string;
using std::vector;

namespace sbsearch
{
    namespace cli
    {
        void conflicting_options(const boost::program_options::variables_map &vm,
                                 const char *opt1, const char *opt2)
        {
            if (vm.count(opt1) && !vm[opt1].defaulted() && vm.count(opt2) && !vm[opt2].defaulted())
                throw std::logic_error(string("Conflicting options '") + opt1 + "' and '" + opt2 + "'.");
        }

        void action_conflicting_option(const boost::program_options::variables_map &vm,
                                       const char *action, const char *option)
        {
            if (!vm.count("action"))
                return;

            if ((vm["action"].as<string>() == action) & vm.count(option))
                throw std::logic_error(string("Action '") + action + "' does not use option '" + option + "'.");
        }

        void option_dependency(const boost::program_options::variables_map &vm,
                               const char *for_what, const char *required_option)
        {
            if (vm.count(for_what) && !vm[for_what].defaulted())
                if (vm.count(required_option) == 0 || vm[required_option].defaulted())
                    throw std::logic_error(string("Option '") + for_what + "' requires option '" + required_option + "'.");
        }

        void action_dependency(const boost::program_options::variables_map &vm,
                               const char *action, const char *required_option)
        {
            if (!vm.count("action"))
                return;

            if ((vm["action"].as<string>() == action) & !vm.count(required_option))
                throw std::logic_error(string("Action '") + action + "' requires option '" + required_option + "'.");
        }

        void validate(boost::any &v,
                      const vector<string> &values,
                      Date *, int)
        {
            using namespace boost::program_options;

            // Make sure no previous assignment was made.
            validators::check_first_occurrence(v);

            // Extract the first string from 'values'. If there is more than
            // one string, it's an error, and exception will be thrown.
            const string &s = validators::get_single_string(values);
            if (s.size() != 10)
                throw validation_error(validation_error::invalid_option_value);

            // Pattern match and convert the date to MJD (double).
            int y, m, d;
            y = (int)std::stoul(s.substr(0, 4).c_str());
            m = (int)std::stoul(s.substr(5, 2).c_str());
            d = (int)std::stoul(s.substr(8, 2).c_str());

            double djm0, djm;
            int status = iauCal2jd(y, m, d, &djm0, &djm);
            if (status)
                throw validation_error(validation_error::invalid_option_value);

            v = boost::any(Date{s, djm});
        }

        bool confirm(const string prompt)
        {
            string response;
            std::cout << prompt << " " << std::flush;
            std::cin >> response;
            return ((response[0] == 'y') | (response[1] == 'Y'));
        }

        const string read_file(const string &file)
        {
            std::ifstream inf(file);
            if (!inf.is_open())
                throw std::runtime_error("failed to open " + file);

            string table;

            inf.seekg(0, std::ios::end);
            table.reserve(inf.tellg());
            inf.seekg(0, std::ios::beg);

            table.assign((std::istreambuf_iterator<char>(inf)),
                         std::istreambuf_iterator<char>());
            return table;
        }

        size_t save_http_data(void *buffer, size_t size, size_t nmemb, void *data)
        {
            string *table = (string *)data;
            size_t realsize = size * nmemb;

            try
            {
                if (table->capacity() < (table->size() + realsize))
                    table->reserve(table->capacity() + realsize);
            }
            catch (std::exception &e)
            {
                cerr << "Cannot reserve string capacity " << table->size() + realsize << ", capacity is already " << table->capacity() << "\n";
                throw;
            }
            table->append((char *)buffer, realsize);
            return realsize;
        }

        const string format_horizons_command(const string target)
        {
            // Temporary comet designation?  C/2001 Q4, P/2003 CC22
            bool comet = (target.find_first_of("CPDI") == 0) & (target[1] == '/');

            // Permanent comet designation? 2P, 2I, 100P.  Check that there is a
            // parseable number without whitespace before the letter.
            if (!comet)
            {
                size_t letter = target.find_first_of("PI");
                if ((letter != string::npos) & !std::isspace(target[letter - 1]))
                {
                    try
                    {
                        std::stoul(target.substr(0, letter));
                        comet = true;
                    }
                    catch (std::exception &e)
                    {
                        // the string could not be converted to a number
                    }
                }
            }

            if (comet)
                return "DES=" + target + ";NOFRAG;CAP;";

            if (target.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ") == string::npos)
                // probably a permanent asteroid ID, just append a semicolon
                return target + ";";

            // otherwise, assume it is a temporary ID
            return "DES=" + target + ";";
        }

        const string from_horizons(const string target, const string center, const Date start_date, const Date stop_date, const string time_step, const bool verbose)
        {
            char query[1024];
            sprintf(query, R"(
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
                    format_horizons_command(target).c_str(),
                    center.c_str(),
                    start_date.ymd.c_str(),
                    stop_date.ymd.c_str(),
                    time_step.c_str());

            if (verbose)
                cerr << query << "\n";

            string table;
            try
            {
                char error_message[CURL_ERROR_SIZE];

                CURL *handle = curl_easy_init();
                curl_easy_setopt(handle, CURLOPT_URL, "https://ssd.jpl.nasa.gov/api/horizons_file.api");
                curl_easy_setopt(handle, CURLOPT_WRITEFUNCTION, save_http_data);
                curl_easy_setopt(handle, CURLOPT_WRITEDATA, (void *)&table);
                curl_easy_setopt(handle, CURLOPT_ERRORBUFFER, error_message);

                curl_mime *multipart = curl_mime_init(handle);
                curl_mimepart *part = curl_mime_addpart(multipart);
                curl_mime_name(part, "format");
                curl_mime_data(part, "text", CURL_ZERO_TERMINATED);
                part = curl_mime_addpart(multipart);
                curl_mime_name(part, "input");
                curl_mime_data(part, query, CURL_ZERO_TERMINATED);
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
                if (verbose)
                    cerr << table;
                throw;
            }

            return table;
        }

        Ephemeris::Data parse_horizons(const string &table)
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
                int n = table.find("\n", row_start) - row_start;
                string line = table.substr(row_start, n);
                row_start += n + 1;

                if ((line == "$$EOE") | (row_start >= data_end))
                    break;

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
            }

            return data;
        }
    }
}