#include <iostream>
#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "cli.h"
#include "config.h"
#include "constants.h"
#include "sbs-configure/arguments.h"

using namespace sbsearch;
using namespace sbsearch::cli;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace sbsearch::sbs_configure
{
    Arguments get_arguments(int argc, char *argv[], Indexer::Options current_options)
    {
        namespace po = boost::program_options;

        Arguments args;
        args.indexer_options = current_options;
        args.reconfigured = false;

        po::options_description options("Options");
        options.add_options()(
            "create,c", po::bool_switch(&args.create), "create database tables and indices if any do not exist")(
            "drop,d", po::bool_switch(&args.drop_indices), "drop observations table indices")(
            "reindex,r", po::bool_switch(&args.reindex), "rebuild the observations spatial indices")(
            "max-spatial-index-cells", po::value<int>(), "maximum number of spatial index cells per observation")(
            "min-spatial-resolution", po::value<double>(), "set minimum spatial level to this angular scale, arcmin")(
            "max-spatial-resolution", po::value<double>(), "set maximum spatial level to this angular scale, arcmin")(
            "min-spatial-level", po::value<int>(), "minimum spatial level")(
            "max-spatial-level", po::value<int>(), "maximum spatial level")(
            "threads", po::value<unsigned int>(&args.threads)->default_value(2), "number of threads for re-indexing");

        po::options_description general = get_common_options((CommonArguments *)&args);

        po::options_description all("");
        all.add(options).add(general);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(all).run(), vm);
        po::notify(vm);

        if (vm.count("version"))
        {
            cout << "SBSearch " << SBSEARCH_VERSION << "\n";
            exit(0);
        }

        if (vm.count("help"))
        {
            cout << "Usage: sbs-configure [options...]\n\n"
                 << "Configure and re-index an sbsearch database.\n\n"
                 << all << "\n";
            exit(0);
        }

        validate_common_options(vm);
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

        args.indexer_options.verify();

        if ((args.threads < 1) || (args.threads > MAX_THREADS))
            throw std::range_error("Number of threads must be between 1 and " + std::to_string(MAX_THREADS) + ".");

        return args;
    }
}