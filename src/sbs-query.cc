#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/json.hpp>
#include <curl/curl.h>
#include <s2/s2latlng.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>

#include "config.h"
#include "date.h"
#include "files.h"
#include "horizons.h"
#include "ephemeris.h"
#include "logging.h"
#include "moving_target.h"
#include "sbsearch.h"
#include "sbsdb/postgresql.h"
#include "cli.h"

using namespace sbsearch;
using namespace sbsearch::cli;
using std::cerr;
using std::cout;
using std::string;
using std::vector;

struct Arguments : CommonArguments
{
    string target;
    bool fixed_target;
    bool input_file;

    vector<string> sources;
    bool parallax;
    bool use_uncertainty;
    double padding = 0;
    double arc_length = 4;
    double time_period = 30;
    bool approximate;
    bool save;
    string info_file;
    string output_filename;
    OutputFormat output_format;
    bool show_fov = false;

    IntersectionType intersection_type = IntersectsArea;

    string file;
    bool horizons;
    bool small_body;

    string observer;
    optional<Date> start_date, stop_date;
    string time_step;

    bool cache;
};

Arguments get_arguments(int argc, char *argv[])
{
    using namespace boost::program_options;

    Arguments args;

    positional_options_description positional;
    positional.add("target", 1);

    options_description hidden("Hidden options");
    hidden.add_options()("target", value<string>(&args.target), "target");

    options_description common_options("Moving / fixed target common options");
    common_options.add_options()(
        "input,i", bool_switch(&args.input_file), "read target names from an input file")(
        "fixed", bool_switch(&args.fixed_target), "indicates <target> is an RA, Dec pair in degrees, e.g., \"123.45 67.890\"")(
        "source,s", value<vector<string>>(&args.sources), "only search this source data set, may be specified multiple times")(
        "padding,p", value<double>(&args.padding), "areal search around query, in arcminutes")(
        "arc-length,arc", value<double>(&args.arc_length), "maximum arc length for ephemeris splitting, degrees")(
        "time-period", value<double>(&args.time_period), "maximum time period for ephemeris splitting, days")(
        "approximate,a", bool_switch(&args.approximate), "return approximate results")(
        "output,o", value<string>(&args.output_filename), "save the results to this file")(
        "format,f", value<OutputFormat>(&args.output_format)->default_value(TableFormat), "output file format: table (default) or json")(
        "show-fov", bool_switch(&args.show_fov), "show fields of view in output table")(
        "info", value<string>(&args.info_file), "save query information to this file, JSON format");

    options_description fixed_target_options("Fixed target options");
    fixed_target_options.add_options()(
        "intersection-type", value<IntersectionType>(&args.intersection_type), "areal search type: ContainsArea, IntersectsArea (default), or ContainedByArea");

    options_description moving_target_options("Moving target options");
    moving_target_options.add_options()(
        "major-body", bool_switch(&args.small_body)->default_value(true), "moving target is a major body")(
        "format-help", "display help on file formats and exit")(
        "file", value<string>(&args.file), "read ephemeris from this file (JSON or Horizons format)")(
        "horizons", bool_switch(&args.horizons), "generate ephemeris with JPL/Horizons")(
        "observer", value<string>(&args.observer)->default_value("500@399"), "observer location for Horizons query")(
        "start", value<optional<Date>>(&args.start_date), "start date for query [YYYY-MM-DD or MJD]")(
        "stop,end", value<optional<Date>>(&args.stop_date), "stop date for query [YYYY-MM-DD or MJD]")(
        "step", value<string>(&args.time_step)->default_value("1d"), "time step size and unit for Horizons query")(
        "use-uncertainty,u", bool_switch(&args.use_uncertainty), "areal search around ephemeris position using the ephemeris uncertainty")(
        "no-cache", bool_switch(&args.cache)->default_value(true), "do not use a file cache for Horizons queries")(
        "no-parallax", bool_switch(&args.parallax)->default_value(true), "do not account for moving target parallax between observatory and the Earth's center")(
        "save", bool_switch(&args.save), "save the results to the found object database");

    options_description general = get_common_options((CommonArguments *)&args);

    options_description visible("");
    visible.add(common_options).add(fixed_target_options).add(moving_target_options).add(general);

    options_description all("");
    all.add(visible).add(hidden);

    variables_map vm;
    boost::program_options::store(command_line_parser(argc, argv).options(all).positional(positional).run(), vm);
    boost::program_options::notify(vm);

    if (vm.count("version"))
    {
        cout << "SBSearch version " << SBSEARCH_VERSION << "\n";
        exit(0);
    }

    if (vm.count("help") | !vm.count("target"))
    {
        cout << "Usage: sbs-query <target> [options...]\n\n"
             << "Find observations of a moving or fixed target.\n\n"
             << "<target> is the ephemeris target name / designation or,\n"
             << "with --fixed-target, an RA, Dec pair.  Use --input to\n"
             << "indicate that <target> is a file listing multiple\n"
             << "targets\n"
             << visible << "\n";

        if (!vm.count("target"))
            cout << "\ntarget is a required argument\n";

        exit(0);
    }

    conflicting_options(vm, "file", "horizons");
    conflicting_options(vm, "file", "observer");
    conflicting_options(vm, "file", "fixed-target");
    conflicting_options(vm, "file", "input");
    conflicting_options(vm, "fixed-target", "horizons");
    conflicting_options(vm, "fixed-target", "parallax");
    conflicting_options(vm, "fixed-target", "use-uncertainty");
    conflicting_options(vm, "fixed-target", "observer");
    option_dependency(vm, "horizons", "start");
    option_dependency(vm, "horizons", "stop");

    return args;
}

template <typename DB>
const Observations query_fixed_target(const Arguments &args, const string &coordinates, SBSearch<DB> &sbs)
{
    // convert target coordinates into S2Point
    const int delimiter = coordinates.find_first_of(", ");
    const double ra = std::stod(coordinates.substr(0, delimiter));
    const double dec = std::stod(coordinates.substr(delimiter + 1));
    S2Point point = S2LatLng::FromDegrees(dec, ra).Normalized().ToPoint();

    // default is to search over all time
    const double mjd_start = args.start_date.value_or(Date(0)).mjd();
    const double mjd_stop = args.stop_date.value_or(Date(100000)).mjd();

    // set options and search
    typename SBSearch<DB>::FindOptions find_options = {.mjd_start = mjd_start,
                                                       .mjd_stop = mjd_stop,
                                                       .padding = args.padding,
                                                       .approximate = args.approximate};
    if (args.padding > 0)
        find_options.intersection_type = args.intersection_type;

    Observations observations;
    observations = sbs.find_observations(point, find_options);

    return observations;
}

template <typename DB>
const Founds query_moving_target(const Arguments &args, const string &designation, SBSearch<DB> &sbs)
{
    // set up moving target
    MovingTarget target = sbsdb::get::moving_target(sbs.db(), designation);

    // default is to search over all time
    const double mjd_start = args.start_date.value_or(Date(0)).mjd();
    const double mjd_stop = args.stop_date.value_or(Date(100000)).mjd();

    Ephemeris eph;
    if (!args.file.empty())
    {
        message("Reading ephemeris from file " + args.file);
        eph = Ephemeris(target, Horizons::parse(read_file(args.file)));
    }
    else if (args.horizons)
    {
        message("Fetching ephemeris for " + target.to_string() + " from Horizons.");
        Horizons horizons(target,
                          args.observer,
                          mjd_start,
                          mjd_stop,
                          args.time_step,
                          args.cache);
        eph = Ephemeris(target, horizons.get_ephemeris_data());
    }
    else
    {
        message("Fetching ephemeris for " + target.to_string() + " from database.");

        eph = sbsdb::get::ephemeris(sbs.db(), target, mjd_start, mjd_stop);
        if (eph.num_vertices() == 0)
            throw std::runtime_error("No ephemeris data for target found in database.");
    }

    // set up search options
    eph.mutable_options()->use_uncertainty = args.use_uncertainty;
    typename SBSearch<DB>::FindOptions find_options = {.mjd_start = mjd_start,
                                                       .mjd_stop = mjd_stop,
                                                       .parallax = args.parallax,
                                                       .save = args.save,
                                                       .padding = args.padding,
                                                       .arc_length = args.arc_length,
                                                       .time_period = args.time_period,
                                                       .approximate = args.approximate,
                                                       .save_info = !args.info_file.empty()};

    // search
    Founds founds;
    if (args.sources.empty())
        founds = sbs.find_observations(eph, find_options);
    else
    {
        for (const string &source : args.sources)
        {
            find_options.source = source;
            founds.append(sbs.find_observations(eph, find_options));
        }
    }

    return founds;
}

template <typename DB>
void sbs_query(int argc, char *argv[])
{
    Arguments args = get_arguments(argc, argv);

    // Set log level
    int log_level = INFO;
    if (args.verbose)
        log_level = DEBUG;

    SBSearch<DB> sbs(args.database, {args.log_file, log_level});
    message("SBSearch moving target query tool.\n");

    // setup target name array
    vector<string> targets;
    if (args.input_file)
    {
        std::ifstream input(args.target);
        for (string line; std::getline(input, line);)
            if ((line.size() > 0) & (line[0] != '#'))
                targets.push_back(line);
    }
    else
        targets.push_back(args.target);

    // Set up output stream: file or stdout
    std::ostream *os;
    std::ofstream outf;
    if (args.output_filename.empty())
        os = &cout;
    else
    {
        outf.open(args.output_filename);
        os = &outf;
    }

    // fixed target search
    if (args.fixed_target)
    {
        Observations observations;
        for (string target : targets)
        {
            Observations new_observations = query_fixed_target(args, target, sbs);
            observations.append(new_observations);
        }

        // output
        if (args.output_format == TableFormat)
        {
            if (observations.size() > 0)
                observations.format.show_fov = args.show_fov;
            *os << observations;
        }
        else
        {
            json::array array;
            for (Observation obs : observations)
                array.emplace_back(obs.as_json());

            *os << array;
        }
    }
    else
    // moving target search
    {
        Founds founds;
        for (string target : targets)
            founds.append(query_moving_target(args, target, sbs));

        cout << "\n";

        // output
        if (args.output_format == TableFormat)
        {
            if (founds.size() > 0)
                founds.data[0].observation.format.show_fov = args.show_fov;
            *os << founds;
        }
        else
            *os << founds.as_json();
    }

    *os << "\n";

    if (!args.info_file.empty())
    {
        std::ofstream outf(args.info_file);
        outf << sbs.query_info();
    }
}

int main(int argc, char *argv[])
{
    curl_global_init(CURL_GLOBAL_ALL);

    try
    {
        // identify which flavor of database to use
        string database = get_arguments(argc, argv).database;
        if (database.substr(0, 8) == "postgres")
            sbs_query<sbsdb::Postgresql>(argc, argv);
        else
            throw DatabaseError("Cannot determine database type from string " + database);
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
