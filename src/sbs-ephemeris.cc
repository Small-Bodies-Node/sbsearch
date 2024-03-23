#include <iostream>
#include <stdexcept>
#include <string>
#include <boost/program_options.hpp>
#include <curl/curl.h>

#include "config.h"
#include "date.h"
#include "ephemeris.h"
#include "files.h"
#include "horizons.h"
#include "logging.h"
#include "moving_target.h"
#include "sbsearch.h"
#include "sbs-cli.h"

using namespace sbsearch;
using namespace sbsearch::cli;
using std::cerr;
using std::cout;
using std::string;

struct Arguments
{
    string action;

    string file;

    string target;
    bool small_body;
    string observer;
    Date start_date, stop_date;
    string time_step;

    bool remove_all;

    bool cache;
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
        "horizons", "generate ephemeris with JPL/Horizons")(
        "observer", value<string>(&args.observer)->default_value("500@399"), "observer location for Horizons query");

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
        "no-cache", bool_switch(&args.cache)->default_value(true), "do not use the cache for Horizons queries")(
        "database,D", value<string>(&args.database)->default_value("sbsearch.db"), "SBSearch database name or file")(
        "db-type,T", value<string>(&args.database_type)->default_value("sqlite3"), "database type")(
        "log-file,L", value<string>(&args.log_file)->default_value("sbsearch.log"), "log file name")(
        "help,h", "display this help and exit")(
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

// add ephemeris data from file or horizons
void add(const Arguments &args, SBSearch &sbs)
{
    MovingTarget target = sbs.db()->get_moving_target(args.target, args.small_body);

    cout << "\nAdding ephemeris for " << target.designation() << ".\n";
    string table;
    Ephemeris::Data data;

    if (!args.file.empty())
    {
        cout << "Reading ephemeris from file " << args.file << ".\n";
        table = read_file(args.file);
        data = Horizons::parse(table);
    }
    else
    {
        cout << "Fetching ephemeris from Horizons API.\n";
        Horizons horizons(
            target,
            args.observer,
            args.start_date,
            args.stop_date,
            args.time_step,
            args.cache);
        data = horizons.get_ephemeris_data();
    }

    Ephemeris eph = Ephemeris(target, data);
    if (eph.num_vertices() == 0)
    {
        cerr << table;
        throw std::runtime_error("Empty ephemeris data.");
    }
    cout << "Read " << eph.num_vertices() << " ephemeris epochs.\n\n";

    sbs.add_ephemeris(eph);
}

// remove the ephemeris points by target, optionally for a date range
void remove(const Arguments &args, SBSearch &sbs)
{
    MovingTarget target = sbs.db()->get_moving_target(args.target);
    if (args.remove_all)
        sbs.db()->remove_ephemeris(target);
    else
        sbs.db()->remove_ephemeris(target, args.start_date.mjd(), args.stop_date.mjd());
}

int main(int argc, char *argv[])
{
    curl_global_init(CURL_GLOBAL_ALL);

    try
    {
        Arguments args = get_arguments(argc, argv);

        // Set log level
        int log_level = sbsearch::INFO;
        if (args.verbose)
            log_level = sbsearch::DEBUG;

        SBSearch sbs(SBSearch::sqlite3, args.database, {args.log_file, log_level});
        Logger::info() << "SBSearch ephemeris management tool." << std::endl;

        if (args.action == "add") // add data to database
            add(args, sbs);
        else if (args.action == "remove") // remove data from database
            remove(args, sbs);
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
