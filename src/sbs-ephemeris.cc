#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <streambuf>
#include <string>
#include <vector>
#include <boost/program_options.hpp>
#include <curl/curl.h>

#include "config.h"
#include "ephemeris.h"
#include "logging.h"
#include "moving_target.h"
#include "sbsearch.h"
#include "sbs-cli.h"
#include "util.h"
#include "sofa/sofa.h"

using sbsearch::Ephemeris;
using sbsearch::Logger;
using sbsearch::MovingTarget;
using sbsearch::SBSearch;
using std::cerr;
using std::cout;
using std::string;
using std::vector;

// Read file contents into a string.
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

// Write HTTP data to a string.
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

// Format the target name as a Horizons query COMMAND string. If the target
// appears to be a comet, then fragment searching is disabled (NOFRAG) and the
// closest apparition is requested (CAP).  Possible temporary asteroidal
// designations (e.g., 2000 XY) are prefixed with "DES="".
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

// Get an ephemeris table  Horizons's API.
const string from_horizons(const string target, const Date start_date, const Date stop_date, const string time_step, const bool verbose)
{
    char query[1024];
    sprintf(query, R"(
!$$SOF
MAKE_EPHEM=YES
COMMAND='%s'
EPHEM_TYPE=OBSERVER
CENTER='500@399'
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

// Parse a Horizons ephemeris table.
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

struct Arguments
{
    string target;
    string action;
    string file;
    Date start_date, stop_date;
    string time_step;
    bool remove_all;
    string database;
    string database_type;
    string log_file;
    bool verbose;
};

Arguments get_arguments(int argc, char *argv[])
{
    using namespace boost::program_options;

    Arguments args;

    positional_options_description positional;
    positional.add("action", 1).add("target", 1);

    options_description hidden("Hidden options");
    hidden.add_options()("action", value<string>(&args.action), "ephemeris action");
    hidden.add_options()("target", value<string>(&args.target), "target name");

    options_description add_options("Options for add action");
    add_options.add_options()(
        "file", value<string>(&args.file), "read ephemeris from this file (JSON or Horizons format)")(
        "horizons", "generate ephemeris with JPL/Horizons");

    options_description remove_options("Options for remove action");
    remove_options.add_options()(
        "all", bool_switch(&args.remove_all), "remove all ephemeris data");

    options_description date_range("Options for date ranges");
    date_range.add_options()(
        "start", value<Date>(&args.start_date),
        "start date for adding/removing ephemeris data [YYYY-MM-DD]")(
        "stop,end", value<Date>(&args.stop_date),
        "stop date for adding/removing ephemeris data [YYYY-MM-DD]")(
        "step", value<string>(&args.time_step)->default_value("1d"), "time step size and unit for Horizons ephemeris generation");

    options_description general("General options");
    general.add_options()(
        "database,D", value<string>(&args.database)->default_value("sbsearch.db"), "SBSearch database name or file")(
        "db-type,T", value<string>(&args.database_type)->default_value("sqlite3"), "database type")(
        "log-file,L", value<string>(&args.log_file)->default_value("sbsearch.log"), "log file name")(
        "help", "display this help and exit")(
        "version", "output version information and exit")(
        "verbose,v", bool_switch(&args.verbose), "show debugging messages");

    options_description visible("");
    visible.add(add_options).add(remove_options).add(date_range).add(general);

    options_description all("");
    all.add(visible).add(hidden);

    options_description add_action("");
    add_action.add(add_options).add(date_range).add(general);

    options_description remove_action("");
    remove_action.add(remove_options).add(date_range).add(general);

    variables_map vm;
    boost::program_options::store(command_line_parser(argc, argv).options(all).positional(positional).run(), vm);
    boost::program_options::notify(vm);

    if (vm.count("version"))
    {
        cout << "SBSearch version " << SBSEARCH_VERSION << "\n";
        exit(0);
    }

    if (vm.count("help") | !vm.count("target") | !vm.count("action"))
    {
        // help for a specific action?
        if (args.action == "add")
        {
            cout << "Usage: sbs-ephemeris add <target> [options...]\n"
                 << "Add ephemeris data for a target to the database.\n\n"
                 << "<target> is the ephemeris target name / designation\n"
                 << add_action << "\n";
        }
        else if (args.action == "remove")
        {
            cout << "Usage: sbs-ephemeris remove <target> [options...]\n"
                 << "Remove ephemeris data for a target from the database.\n\n"
                 << "<target> is the ephemeris target name / designation\n"
                 << remove_action << "\n";
        }
        else
        {
            cout << "Usage: sbs-ephemeris <action> <target> [options...]\n\n"
                 << "Manage sbsearch ephemeris data.\n\n"
                 << "<action> is one of {add, remove}\n"
                 << "<target> is the ephemeris target name / designation\n"
                 << visible << "\n";
        }

        if ((args.action == "add") | (args.action.empty()))
        {
            cout << "Horizons ephemeris files require the CSV format, angles formatted in degrees,\n"
                 << "dates as Julian days, range in au, use the ICRF reference frame, and\n"
                 << "observer quantities 1, 9, 19, 20, 23, 24, 27, 37, and 41.  Extra precision\n"
                 << "is optional.\n";
        }

        if (!vm.count("action") | !vm.count("target"))
            cout << "\naction and target are required arguments\n";

        exit(0);
    }

    conflicting_options(vm, "file", "horizons");
    option_dependency(vm, "horizons", "start");
    option_dependency(vm, "start", "stop");

    if ((args.action == "remove") & (!args.remove_all) & (!vm.count("start")))
        throw std::logic_error("remove action requires a date range or --all");

    return args;
}

int main(int argc, char *argv[])
{
    curl_global_init(CURL_GLOBAL_ALL);

    try
    {
        Arguments args = get_arguments(argc, argv);
        cout << "SBSearch ephemeris management tool.\n\n";

        SBSearch sbs(SBSearch::sqlite3, args.database, args.log_file);

        // Set log level
        if (args.verbose)
            Logger::get_logger().log_level(sbsearch::DEBUG);
        else
            Logger::get_logger().log_level(sbsearch::INFO);

        if (args.action == "add") // add data to database
        {
            MovingTarget target = sbs.get_moving_target(args.target);

            cout << "\nAdding ephemeris for " << target.designation() << ".\n";
            string table;

            if (!args.file.empty())
            {
                cout << "Reading ephemeris from file " << args.file << ".\n";
                table = read_file(args.file);
            }
            else
            {
                cout << "Fetching ephemeris from Horizons API.\n";
                table = from_horizons(args.target, args.start_date, args.stop_date, args.time_step, args.verbose);
            }

            Ephemeris eph = Ephemeris{MovingTarget(args.target), parse_horizons(table)};
            if (eph.num_vertices() == 0)
            {
                cerr << table;
                throw std::runtime_error("Empty ephemeris data.");
            }
            cout << "Read " << eph.num_vertices() << " ephemeris epochs.\n\n";

            sbs.add_ephemeris(eph);
        }
        else if (args.action == "remove") // remove data from database
        {
            MovingTarget target = sbs.get_moving_target(args.target);
            if (args.remove_all)
                sbs.remove_ephemeris(target);
            else
                sbs.remove_ephemeris(target, args.start_date.mjd, args.stop_date.mjd);
        }
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
