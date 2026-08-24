#include <iostream>
#include <string>
#include <vector>

#include "config.h"
#include "constants.h"
#include "env.h"
#include "exceptions.h"
#include "indexer.h"
#include "logging.h"
#include "sbsearch.h"
#include "sbs-configure/arguments.h"
#include "sbsdb/postgresql.h"
#include "cli.h"

using namespace sbsearch;
using namespace sbsearch::cli;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace sbsearch::sbs_configure
{
    template <typename DB>
    void sbs_configure(int argc, char **argv)
    {
        Arguments args = get_arguments(argc, argv);

        SBSearch<DB> sbs(args.database, {args.log_file,
                                         args.log_level(),
                                         args.create,
                                         static_cast<unsigned char>(args.threads)});
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
             << "\n\n";

        if (args.drop_indices)
        {
            cout << "Dropping observations indices." << endl;
            sbs.drop_observations_indices();
        }

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
                 << "\n\n";
        }

        if (args.reindex || args.reconfigured)
        {
            cout << "Re-indexing..." << endl;
            sbs.reindex_database_terms(args.indexer_options);
        }
    }
}

int main(int argc, char **argv)
{
    try
    {
        // get basic CLI stuff first to tell us which flavor of database to use
        string database = sbs_configure::get_arguments(argc, argv).database;
        if (database.substr(0, 8) == "postgres")
            sbs_configure::sbs_configure<sbsdb::Postgresql>(argc, argv);
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