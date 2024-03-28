#include <iostream>
#include <string>
#include <vector>
#include <boost/program_options.hpp>
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
#include "sbs-cli.h"
#include "util.h"

using namespace sbsearch;
using namespace sbsearch::cli;
using std::cerr;
using std::cout;
using std::string;
using std::vector;

struct Arguments
{
    string target;
    bool fixed_target;
    bool input_file;

    vector<string> sources;
    bool parallax;
    bool use_uncertainty;
    double padding = 0;
    bool save;
    string output_filename;

    string file;
    bool horizons;
    bool small_body;

    string observer;
    Date start_date, stop_date;
    string time_step;

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
    positional.add("target", 1);

    options_description hidden("Hidden options");
    hidden.add_options()("target", value<string>(&args.target), "target");

    options_description common_options("Moving / fixed target common options");
    common_options.add_options()(
        "input,i", bool_switch(&args.input_file), "read target names from an input file")(
        "fixed-target,f", bool_switch(&args.fixed_target), "indicates <target> is an RA, Dec pair in degrees, e.g., \"123.45 67.890\"")(
        "source,s", value<vector<string>>(&args.sources), "only search this source data set, may be specified multiple times")(
        "padding,p", value<double>(&args.padding), "areal search around query, in arcsec")(
        "output,o", value<string>(&args.output_filename), "save the results to this file");

    options_description moving_target_options("Moving target options");
    moving_target_options.add_options()(
        "major-body", bool_switch(&args.small_body)->default_value(true), "moving target is a major body")(
        "format-help", "display help on file formats and exit")(
        "file", value<string>(&args.file), "read ephemeris from this file (JSON or Horizons format)")(
        "horizons", bool_switch(&args.horizons), "generate ephemeris with JPL/Horizons")(
        "observer", value<string>(&args.observer)->default_value("500@399"), "observer location for Horizons query")(
        "start", value<Date>(&args.start_date), "start date for query [YYYY-MM-DD]")(
        "stop,end", value<Date>(&args.stop_date), "stop date for query [YYYY-MM-DD]")(
        "step", value<string>(&args.time_step)->default_value("1d"), "time step size and unit for Horizons query")(
        "use-uncertainty,u", bool_switch(&args.use_uncertainty), "areal search around ephemeris position using the ephemeris uncertainty")(
        "no-cache", bool_switch(&args.cache)->default_value(true), "do not use a file cache for Horizons queries")(
        "no-parallax", bool_switch(&args.parallax)->default_value(true), "do not account for moving target parallax between observatory and the Earth's center")(
        "save", bool_switch(&args.save), "save the results to the found object database");

    options_description general("General options");
    general.add_options()(
        "database,D", value<string>(&args.database)->default_value("sbsearch.db"), "SBSearch database name or file")(
        "db-type,T", value<string>(&args.database_type)->default_value("sqlite3"), "database type")(
        "log-file,L", value<string>(&args.log_file)->default_value("sbsearch.log"), "log file name")(
        "help,h", "display this help and exit")(
        "version", "output version information and exit")(
        "verbose,v", bool_switch(&args.verbose), "show debugging messages");

    options_description visible("");
    visible.add(common_options).add(moving_target_options).add(general);

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

void query_fixed_target(const Arguments &args, SBSearch &sbs)
{
    const int delimiter = args.target.find_first_of(", ");
    const double ra = std::stod(args.target.substr(0, delimiter));
    const double dec = std::stod(args.target.substr(delimiter + 1));
    S2LatLng latlng = S2LatLng::FromDegrees(dec, ra).Normalized();

    // default is to search everything, but the user may limit the query
    double mjd_start = (args.start_date.mjd() == UNDEF_TIME) ? 0 : args.start_date.mjd();
    double mjd_stop = (args.stop_date.mjd() == UNDEF_TIME) ? 100000 : args.stop_date.mjd();

    SBSearch::SearchOptions options = {.mjd_start = mjd_start,
                                       .mjd_stop = mjd_stop};

    Observations obs;
    if (args.padding > 0)
    {
        double r = args.padding * ARCSEC;
        cerr << "here\n";
        vector<S2LatLng> latlngs = ellipse(32, latlng, r, r, 0);
        cerr << "here\n";

        vector<S2Point> points(latlngs.size());
        std::transform(latlngs.begin(), latlngs.end(), points.begin(),
                       [](const S2LatLng &ll)
                       { return ll.ToPoint(); });
        cerr << "here\n";

        S2Polygon polygon;
        makePolygon(points, polygon);

        obs = sbs.find_observations(polygon, options);
    }
    else
        obs = sbs.find_observations(latlng.ToPoint(), options);

    cout << obs;
}

void query_moving_target(const Arguments &args, SBSearch &sbs)
{
    MovingTarget target = sbs.db()->get_moving_target(args.target);
    Ephemeris eph;

    // default is to search everything, but the user may limit the query
    double mjd_start = (args.start_date.mjd() == UNDEF_TIME) ? 0 : args.start_date.mjd();
    double mjd_stop = (args.stop_date.mjd() == UNDEF_TIME) ? 100000 : args.stop_date.mjd();

    cout << "\n";
    if (!args.file.empty())
    {
        cout << "Reading ephemeris from file " << args.file << "\n";
        eph = Ephemeris(target, Horizons::parse(read_file(args.file)));
    }
    else if (args.horizons)
    {
        Logger::info() << "Fetching ephemeris for " << target << " from Horizons." << std::endl;
        Horizons horizons(target,
                          args.observer,
                          args.start_date,
                          args.stop_date,
                          args.time_step,
                          args.cache);
        eph = Ephemeris(target, horizons.get_ephemeris_data());
    }
    else
    {
        Logger::info() << "Fetching ephemeris for " << target << " from database." << std::endl;

        eph = sbs.db()->get_ephemeris(target, mjd_start, mjd_stop);
        if (eph.num_vertices() == 0)
            throw std::runtime_error("No ephemeris data for target found in database.");
    }
    cout << "\n";

    eph.mutable_options()->use_uncertainty = args.use_uncertainty;
    eph.mutable_options()->padding = args.padding;
    SBSearch::SearchOptions search_options = {.mjd_start = mjd_start,
                                              .mjd_stop = mjd_stop,
                                              .parallax = args.parallax,
                                              .save = args.save};

    Founds founds;
    if (args.sources.empty())
        founds = sbs.find_observations(eph, search_options);
    else
    {
        for (const string &source : args.sources)
        {
            search_options.source = source;
            founds.append(sbs.find_observations(eph, search_options));
        }
    }
    cout << founds;
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
        Logger::info() << "SBSearch moving target query tool." << std::endl;

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

        for (string target : targets)
        {
            args.target = target;
            if (args.fixed_target)
                query_fixed_target(args, sbs);
            else
                query_moving_target(args, sbs);
        }
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
