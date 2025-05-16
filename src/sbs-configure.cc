// Configure sbsearch database
#include "config.h"

#include <iostream>
#include <string>
#include <vector>

#include "constants.h"
#include "env.h"
#include "exceptions.h"
#include "indexer.h"
#include "logging.h"
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
    bool create;
    Indexer::Options indexer_options;
    bool reconfigured;
    bool reindex;
};

Arguments get_arguments(int argc, char *argv[], Indexer::Options current_options = Indexer::Options())
{
    namespace po = boost::program_options;

    Arguments args;
    args.indexer_options = current_options;
    args.reconfigured = false;

    po::options_description options("Options");
    options.add_options()(
        "create,c", po::bool_switch(&args.create), "create database tables if any do not exist")(
        "reindex,r", po::bool_switch(&args.reindex), "reindex the observations table")(
        "max-spatial-index-cells", po::value<int>(), "maximum number of spatial index cells per observation")(
        "min-spatial-resolution", po::value<double>(), "set minimum spatial level to this angular scale, arcmin")(
        "max-spatial-resolution", po::value<double>(), "set maximum spatial level to this angular scale, arcmin")(
        "min-spatial-level", po::value<int>(), "minimum spatial level")(
        "max-spatial-level", po::value<int>(), "maximum spatial level")(
        "temporal-resolution", po::value<int>(), "temporal resolution, inverse days");

    po::options_description general = get_common_options((CommonArguments *)&args);

    po::options_description all("");
    all.add(options).add(general);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(all).run(), vm);
    po::notify(vm);

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

    if (vm.count("max-spatial-index-cells"))
    {
        args.indexer_options.max_spatial_index_cells(vm["max-spatial-index-cells"].as<int>());
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

template <typename DB>
void sbs_configure(int argc, char **argv)
{
    Arguments args = get_arguments(argc, argv);

    // Set log level
    int log_level = INFO;
    if (args.verbose)
        log_level = DEBUG;

    SBSearch<DB> sbs(args.database, {args.log_file, log_level, args.create});
    Logger::info() << "SBSearch database configuration tool." << std::endl;

    Indexer::Options previous_options = sbs.indexer_options();

    // now get the indexer options, using the previous values as the default
    args = get_arguments(argc, argv, previous_options);

    cout << "\nCurrent index setup:"
         << "\n  Maximum spatial cells: " << previous_options.max_spatial_index_cells()
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
             << "\n  Maximum spatial cells: " << args.indexer_options.max_spatial_index_cells()
             << "\n  Minimum spatial level: "
             << args.indexer_options.min_spatial_level()
             << " (" << args.indexer_options.max_spatial_resolution() / DEG << " deg)"
             << "\n  Maximum spatial level: "
             << args.indexer_options.max_spatial_level()
             << " (" << args.indexer_options.min_spatial_resolution() / DEG << " deg)"
             << "\n  Temporal resolution (1/day): " << args.indexer_options.temporal_resolution()
             << "\n\n";
    }

    if (args.reindex & !args.reconfigured)
        cout << "Re-indexing";

    if (args.reconfigured | args.reindex)
        sbs.reindex(args.indexer_options);
}

int main(int argc, char **argv)
{
    try
    {
        // get basic CLI stuff first to tell us which flavor of database to use
        string database = get_arguments(argc, argv).database;
        if (database.substr(0, 8) == "postgres")
            sbs_configure<sbsdb::Postgresql>(argc, argv);
        else
            throw DatabaseError("Cannot determine database type from string " + database);
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
    return 0;
}