// Configure sbsearch database
#include "config.h"

#include <iostream>
#include <string>

#include "config.h"
#include "indexer.h"
#include "logging.h"
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
    Indexer::Options indexer_options;
    bool reconfigured;

    string database;
    bool create;
    string database_type;
    string log_file;
    bool verbose;
};

Arguments get_arguments(int argc, char *argv[])
{
    using namespace boost::program_options;

    Arguments args;
    args.reconfigured = false;

    options_description options("Options");
    options.add_options()(
        "max-spatial-cells", value<int>(), "maximum number of spatial cells per observation")(
        "min-spatial-resolution", value<double>(), "set minimum spatial level to this angular scale, arcmin")(
        "max-spatial-resolution", value<double>(), "set maximum spatial level to this angular scale, arcmin")(
        "min-spatial-level", value<int>(), "minimum spatial level")(
        "max-spatial-level", value<int>(), "maximum spatial level")(
        "temporal-resolution", value<int>(), "temporal resolution, inverse days");

    options_description general("General options");
    general.add_options()(
        "database,D", value<string>(&args.database)->default_value("sbsearch.db"), "SBSearch database name or file")(
        "create,c", bool_switch(&args.create), "create database if it does not exist")(
        "db-type,T", value<string>(&args.database_type)->default_value("sqlite3"), "database type")(
        "log-file,L", value<string>(&args.log_file)->default_value("sbsearch.log"), "log file name")(
        "help,h", "display this help and exit")(
        "version", "output version information and exit")(
        "verbose,v", bool_switch(&args.verbose), "show debugging messages");

    options_description all("");
    all.add(options).add(general);

    variables_map vm;
    boost::program_options::store(command_line_parser(argc, argv).options(all).run(), vm);
    boost::program_options::notify(vm);

    if (vm.count("version"))
    {
        cout << "SBSearch version " << SBSEARCH_VERSION << "\n";
        exit(0);
    }

    if (vm.count("help"))
    {
        cout << "Usage: sbs-configure [options...]\n\n"
             << "Configure and re-index an sbsearch database.\n\n"
             << all << "\n";
        exit(0);
    }

    conflicting_options(vm, "min-spatial-resolution", "min-spatial-level");
    conflicting_options(vm, "max-spatial-resolution", "max-spatial-level");

    if (vm.count("max-spatial-cells"))
    {
        args.indexer_options.max_spatial_cells(vm["max-spatial-cells"].as<int>());
        args.reconfigured = true;
    }
    if (vm.count("min-spatial-resolution"))
    {
        args.indexer_options.min_spatial_resolution(vm["min-spatial-resolution"].as<double>() * ARCMIN);
        args.reconfigured = true;
    }
    if (vm.count("max-spatial-resolution"))
    {
        args.indexer_options.max_spatial_resolution(vm["max-spatial-resolution"].as<double>() * ARCMIN);
        args.reconfigured = true;
    }
    if (vm.count("min-spatial-level"))
    {
        args.indexer_options.min_spatial_level(vm["min-spatial-level"].as<int>());
        args.reconfigured = true;
    }
    if (vm.count("max-spatial-level"))
    {
        args.indexer_options.max_spatial_level(vm["max-spatial-level"].as<int>());
        args.reconfigured = true;
    }
    if (vm.count("temporal-resolution"))
    {
        args.indexer_options.temporal_resolution(vm["temporal-resolution"].as<int>());
        args.reconfigured = true;
    }

    return args;
}

int main(int argc, char **argv)
{
    try
    {
        Arguments args = get_arguments(argc, argv);

        // Set log level
        int log_level = INFO;
        if (args.verbose)
            log_level = DEBUG;

        SBSearch sbs(SBSearch::sqlite3, args.database, {args.log_file, log_level, args.create});
        Logger::info() << "SBSearch database configuration tool." << std::endl;

        Indexer::Options previous_options = sbs.indexer_options();

        cout << "\nCurrent index setup:"
             << "\n  Maximum spatial cells: " << previous_options.max_spatial_cells()
             << "\n  Minimum spatial level: "
             << previous_options.min_spatial_level()
             << " (" << previous_options.max_spatial_resolution() / DEG << " deg)"
             << "\n  Maximum spatial level: "
             << previous_options.max_spatial_level()
             << " (" << previous_options.min_spatial_resolution() / DEG << " deg)"
             << "\n  Temporal resolution (1/day): " << previous_options.temporal_resolution()
             << "\n\n";
        if (args.reconfigured)
        {
            cout << "\nNew index setup:"
                 << "\n  Maximum spatial cells: " << args.indexer_options.max_spatial_cells()
                 << "\n  Maximum spatial resolution (deg) / level: "
                 << args.indexer_options.max_spatial_resolution() / DEG
                 << " / " << args.indexer_options.min_spatial_level()
                 << "\n  Minimum spatial resolution (deg) / level: "
                 << args.indexer_options.min_spatial_resolution() / DEG
                 << " / " << args.indexer_options.max_spatial_level()
                 << "\n  Temporal resolution (1/day): " << args.indexer_options.temporal_resolution()
                 << "\n\n";
            sbs.reindex(args.indexer_options);
        }
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
    return 0;
}