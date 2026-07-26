#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <boost/program_options.hpp>
#include <curl/curl.h>

#include "config.h"
#include "date.h"
#include "ephemeris.h"
#include "files.h"
#include "horizons.h"
#include "logging.h"
#include "moving_target.h"
#include "orbital_elements.h"
#include "sbsearch.h"
#include "sbsdb.h"
#include "cli.h"

using namespace sbsearch;
using namespace sbsearch::cli;
using std::cerr;
using std::cout;
using std::endl;
using std::string;

struct Arguments : CommonArguments
{
    string action;

    string file;
    bool input_file;

    string target;
    bool major_body;
    string observer;
    optional<Date> start_date, stop_date;
    string time_step;

    double interpolate = -1;
    string output_filename;
    OutputFormat output_format = TABLE;
    DateFormat date_format = DateFormat::MJD;

    bool remove_all;
    bool cache;
};

Arguments get_arguments(int argc, char *argv[])
{
    using namespace boost::program_options;

    Arguments args;

    positional_options_description positional;
    positional.add("action", 1).add("target", 1);

    options_description hidden("Hidden options");
    hidden.add_options()("action", value<string>(&args.action), "ephemeris action");
    hidden.add_options()("target", value<string>(&args.target), "target");

    options_description target_options("Target options");
    target_options.add_options()(
        "input,i", bool_switch(&args.input_file), "read target names from an input file");

    options_description date_range("Options for date ranges");
    date_range.add_options()(
        "start", value<optional<Date>>(&args.start_date),
        "start date for adding/removing ephemeris data [YYYY-MM-DD or MJD]")(
        "stop,end", value<optional<Date>>(&args.stop_date),
        "stop date for adding/removing ephemeris data [YYYY-MM-DD or MJD]")(
        "step", value<string>(&args.time_step)->default_value("VAR 360"), "time step for Horizons ephemeris generation: length and time unit, or \"VAR X\" where X is an angular distance in arcsec");

    options_description list_options("Options for list action");
    list_options.add_options()(
        "interpolate", value<double>(&args.interpolate), "interpolate the database ephemeris to this time step in days")(
        "output,o", value<string>(&args.output_filename), "save ephemeris to this file")(
        "format,f", value<OutputFormat>(&args.output_format), "output file format: table (default) or json")(
        "date", value<DateFormat>(&args.date_format), "date format: mjd (default) or calendar");

    options_description add_options("Options for add action");
    add_options.add_options()(
        "file", value<string>(&args.file), "read ephemeris from this file (JSON or Horizons format)")(
        "horizons", "generate ephemeris with JPL/Horizons")(
        "observer", value<string>(&args.observer)->default_value("500@399"), "observer location for Horizons query")(
        "major-body", bool_switch(&args.major_body), "moving target is a major body");

    options_description remove_options("Options for remove action");
    remove_options.add_options()(
        "all", bool_switch(&args.remove_all), "remove all ephemeris data");

    options_description general = get_common_options((CommonArguments *)&args);

    options_description visible("");
    visible.add(target_options).add(date_range).add(add_options).add(list_options).add(remove_options).add(general);

    options_description all("");
    all.add(visible).add(hidden);

    options_description add_action("");
    add_action.add(add_options).add(date_range).add(general);

    options_description list_action("");
    list_action.add(date_range).add(list_options).add(general);

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

    if (vm.count("help") || !vm.count("target") || !vm.count("action"))
    {
        // help for a specific action?
        if (args.action == "add")
        {
            cout << "Usage: sbs-ephemeris add <target> [options...]\n"
                 << "Add ephemeris data for a target to the database.\n\n"
                 << "<target> is the ephemeris target name / designation\n"
                 << "or, with -i, an input file listing multiple targets\n"
                 << add_action << "\n";
        }
        else if (args.action == "list")
        {
            cout << "Usage: sbs-ephemeris list <target> [options...]\n"
                 << "List ephemeris data for a target in the database.\n\n"
                 << list_action << "\n";
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
                 << "<action> is one of {add, list, remove}\n"
                 << "<target> is the ephemeris target name / designation\n"
                 << visible << "\n";
        }

        if ((args.action == "add") || (args.action.empty()))
        {
            cout << "Horizons ephemeris files require the CSV format, angles formatted in degrees,\n"
                 << "dates as Julian days, range in au, use the ICRF reference frame, and\n"
                 << "observer quantities 1, 9, 19, 20, 23, 24, 27, 37, 41, and 47.  Extra\n"
                 << "precision is optional.\n";
        }

        if (!vm.count("action") || !vm.count("target"))
            cout << "\naction and target are required arguments\n";

        exit(0);
    }

    conflicting_options(vm, "file", "horizons");
    conflicting_options(vm, "file", "input");
    option_dependency(vm, "horizons", "start");
    option_dependency(vm, "start", "stop");

    if ((args.action == "list") && (args.input_file))
        throw std::logic_error("list action and --input-file are not compatible");

    if ((args.action == "remove") && (!args.remove_all) && (!vm.count("start")))
        throw std::logic_error("remove action requires a date range or --all");

    return args;
}

// add ephemeris data from file or horizons
template <typename DB>
void add(const Arguments &args, SBSearch<DB> &sbs)
{
    MovingTarget target = sbsdb::get::moving_target(sbs.db(), args.target, !args.major_body);
    std::cerr << target << std::endl;
    if (!target.moving_target_id())
    {
        sbsdb::add::moving_target(sbs.db(), target);
        cout << "Added moving target " << target.designation()
             << " to the database with ID " << target.moving_target_id().value()
             << "." << endl;
    }

    int count = 0;
    if (!args.file.empty())
    {
        cout << "Reading ephemeris from file " << args.file << ".\n";
        string table = read_file(args.file);
        Ephemeris eph = Ephemeris(target, Horizons::parse(table));
        sbs.add_ephemeris(eph);
        count = eph.num_vertices();
    }
    else
    {
        cout << "Fetching ephemeris from Horizons API for " << target.designation() << " from "
             << args.start_date.value().iso() << " to " << args.stop_date.value().iso()
             << "." << endl;

        // request a maximum of 1 year at a time for VAR steps, 5 years otherwise
        const int years = (args.time_step.find("VAR") == string::npos) ? 5 : 1;
        auto dates = date_ranges(args.start_date.value(), args.stop_date.value(), years * 365.);

        bool first = true;
        Ephemeris ephemeris(target, {});
        for (auto const &[start, stop] : dates)
        {
            Horizons horizons(target,
                              args.observer,
                              start,
                              stop,
                              args.time_step,
                              args.cache);
            Ephemeris eph(target, horizons.get_ephemeris_data());

            // always verify the last item when using VAR steps
            if ((args.time_step.find("VAR") != string::npos) && (eph.data(-1).mjd < stop.mjd()))
            {
                Logger::warning() << "Variable time step detected, but end date was not returned.  Retrying." << endl;
                horizons.time_step("1");
                auto data = horizons.get_ephemeris_data();
                eph.append(Ephemeris(target, data)[1]);
            }

            // trim the first point... it was already added on the last step
            if (!first)
                eph = eph.slice(1);

            ephemeris.append(eph);

            first = false;
        }
        count = ephemeris.num_vertices();
        sbs.add_ephemeris(ephemeris);
    }

    if (count == 0)
        throw std::runtime_error("Empty ephemeris data.");
    else
        cout << "Added " << count << " ephemeris epochs.\n\n";
}

// list ephemeris data in the database, optionally for a date range, optionally
// interpolating, optionally for an observatory
template <typename DB>
void list(const Arguments &args, SBSearch<DB> &sbs)
{
    MovingTarget target = sbsdb::get::moving_target(sbs.db(), args.target, !args.major_body);
    if (!target.moving_target_id())
        throw MovingTargetError("Target not found.");

    Ephemeris eph;
    if (args.interpolate > 0)
    {
        // For interpolation, get the whole ephemeris from the database and then
        // interpolate
        eph = sbsdb::get::ephemeris(sbs.db(), target);

        // never extrapolate
        double start_mjd = std::max(args.start_date.value().mjd(), eph.data(0).mjd.value());
        double stop_mjd = std::min(args.stop_date.value().mjd(), eph.data(-1).mjd.value());

        Ephemeris interpolated;
        interpolated.target(target);
        for (double mjd = start_mjd; mjd < stop_mjd; mjd += args.interpolate)
            interpolated.append(eph.interpolate(mjd));

        interpolated.append(eph.interpolate(stop_mjd));
        eph = interpolated;
    }
    else
    {
        // If we are not interpolating, then search the database based on
        // start/stop dates
        eph = sbsdb::get::ephemeris(
            sbs.db(),
            target,
            args.start_date.value().mjd(),
            args.stop_date.value().mjd());
    }

    if (args.observer != "500@399")
    {
        Observatory observatory = sbsdb::get::observatory(sbs.db(), args.observer);
        eph = eph.parallax_offset(observatory);
    }

    std::ostream *os;
    std::ofstream outf;
    if (args.output_filename.empty())
        os = &cout;
    else
    {
        outf.open(args.output_filename, std::ios::trunc);
        os = &outf;
    }

    // output format
    eph.format.date = args.date_format;
    if (args.output_format == TABLE)
        *os << eph;
    else
    {
        json::array array = eph.as_json();
        *os << array;
    }

    // close file stream
    if (outf.is_open())
        outf.close();
}

// remove the ephemeris points by target, optionally for a date range
template <typename DB>
void remove(const Arguments &args, SBSearch<DB> &sbs)
{
    MovingTarget target = sbsdb::get::moving_target(sbs.db(), args.target, !args.major_body);
    int count;
    if (args.remove_all)
        count = sbsdb::remove::ephemeris(sbs.db(), target);
    else
        count = sbsdb::remove::ephemeris(
            sbs.db(),
            target,
            args.start_date.value().mjd(),
            args.stop_date.value().mjd());

    cout << count << " ephemeris rows removed." << endl;
}

template <typename DB>
void sbs_ephemeris(int argc, char *argv[])
{
    Arguments args = get_arguments(argc, argv);

    SBSearch<DB> sbs(args.database, {args.log_file, args.log_level()});

    Logger::info() << "SBSearch ephemeris management tool." << endl;

    // ephemeris date range is from command line or observations date range
    auto range = sbsdb::get::observations_date_range(sbs.db());
    if ((args.action != "remove") && (!args.start_date || !args.stop_date) && (!range.first || !range.second))
        throw EphemerisError("Observations database is empty: --start and --stop are required.");

    args.start_date = args.start_date ? args.start_date.value() : range.first.value();
    args.stop_date = args.stop_date ? args.stop_date.value() : range.second.value();

    vector<string> targets;
    if (args.input_file)
    {
        std::ifstream input(args.target);
        for (string line; std::getline(input, line);)
            if ((line.size() > 0) && (line[0] != '#'))
                targets.push_back(line);
    }
    else
        targets.push_back(args.target);

    for (string target : targets)
    {
        args.target = target;
        if (args.action == "add") // add data to database
            add(args, sbs);
        else if (args.action == "list") // list ephemeris data
            list(args, sbs);
        else if (args.action == "remove") // remove data from database
            remove(args, sbs);
    }
}

int main(int argc, char *argv[])
{
    curl_global_init(CURL_GLOBAL_ALL);

    try
    {
        // get basic CLI stuff first to tell us which flavor of database to use
        string database = get_arguments(argc, argv).database;
        if (database.substr(0, 8) == "postgres")
            sbs_ephemeris<sbsdb::Postgresql>(argc, argv);
        else
            throw DatabaseError("Cannot determine database type from string " + database);
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
