#include <iostream>
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

using namespace sbsearch;
using namespace sbsearch::cli;
using std::cerr;
using std::cout;
using std::string;
using std::vector;

struct Arguments
{
    string target;
    bool save;
    string output_filename;
    vector<string> sources;
    bool parallax;
    bool use_uncertainty;
    double padding = 0;

    string file;
    bool horizons;

    string observer;
    Date start_date, stop_date;
    string time_step;

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
    hidden.add_options()("target", value<string>(&args.target), "target name");

    options_description options("Options");
    options.add_options()(
        "source,s", value<vector<string>>(&args.sources), "only search this source, may be specified multiple times")(
        "parallax", bool_switch(&args.parallax), "search accounting for parallax between observatory and the Earth's center")(
        "use-uncertainty,u", bool_switch(&args.use_uncertainty), "areal search around ephemeris position using the ephemeris uncertainty")(
        "padding,p", value<double>(&args.padding), "areal search around ephemeris position, in arcsec")(
        "save", bool_switch(&args.save), "save the results to the found object database")(
        "o", value<string>(&args.output_filename), "save the results this file")(
        "file", value<string>(&args.file), "read ephemeris from this file (JSON or Horizons format)")(
        "format-help", "display help on file formats and exit")(
        "horizons", bool_switch(&args.horizons), "generate ephemeris with JPL/Horizons")(
        "observer", value<string>(&args.observer)->default_value("500@399"), "observer location for Horizons query")(
        "start", value<Date>(&args.start_date), "start date for query [YYYY-MM-DD]")(
        "stop,end", value<Date>(&args.stop_date), "stop date for query [YYYY-MM-DD]")(
        "step", value<string>(&args.time_step)->default_value("1d"), "time step size and unit for Horizons query");

    options_description general("General options");
    general.add_options()(
        "database,D", value<string>(&args.database)->default_value("sbsearch.db"), "SBSearch database name or file")(
        "db-type,T", value<string>(&args.database_type)->default_value("sqlite3"), "database type")(
        "log-file,L", value<string>(&args.log_file)->default_value("sbsearch.log"), "log file name")(
        "help,h", "display this help and exit")(
        "version", "output version information and exit")(
        "verbose,v", bool_switch(&args.verbose), "show debugging messages");

    options_description visible("");
    visible.add(options).add(general);

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
        cout << "Usage: sbs-query-moving-target <target> [options...]\n\n"
             << "Find observations of a moving target.\n\n"
             << "<target> is the ephemeris target name / designation\n"
             << visible << "\n";

        if (!vm.count("target"))
            cout << "\ntarget is a required argument\n";

        exit(0);
    }

    conflicting_options(vm, "file", "horizons");
    option_dependency(vm, "horizons", "start");
    option_dependency(vm, "start", "stop");

    return args;
}

void query(const Arguments &args, SBSearch &sbs)
{
    MovingTarget target = sbs.db()->get_moving_target(args.target);
    Ephemeris eph;

    // default is to search everything, but the user may limit the query
    double mjd_start = (args.start_date.mjd == UNDEF_TIME) ? 0 : args.start_date.mjd;
    double mjd_stop = (args.stop_date.mjd == UNDEF_TIME) ? 100000 : args.stop_date.mjd;

    cout << "\n";
    if (!args.file.empty())
    {
        cout << "Reading ephemeris from file " << args.file << "\n";
        eph = Ephemeris(target, parse_horizons(read_file(args.file)));
    }
    else if (args.horizons)
    {
        cout << "Fetching ephemeris for " << target << " from Horizons.\n";
        const string table = from_horizons(target.designation(),
                                           args.observer,
                                           args.start_date,
                                           args.stop_date,
                                           args.time_step,
                                           args.verbose);
        eph = Ephemeris(target, parse_horizons(table));
    }
    else
    {
        cout << "Fetching ephemeris for " << target << " from database.\n";

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

    vector<Found> founds;
    if (args.sources.empty())
        founds = sbs.find_observations(eph, search_options);
    else
    {
        for (const string &source : args.sources)
        {
            search_options.source = source;
            for (const Found &f : sbs.find_observations(eph, search_options))
                founds.push_back(f);
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

        query(args, sbs);
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
