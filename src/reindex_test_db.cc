#include "config.h"

#include <iostream>

#include "indexer.h"
#include "logging.h"
#include "sbsearch.h"
#include "test_db.h"

using sbsearch::Indexer;
using sbsearch::Logger;
using sbsearch::SBSearch;

int main(int argc, char **argv)
{
    SBSearch sbs(SBSearch::sqlite3, "sbsearch_test.db", {.log_file = "sbsearch_test.log"});

    Indexer::Options options;
    options.max_spatial_cells(MAX_SPATIAL_CELLS);
    options.max_spatial_resolution(MAX_SPATIAL_RESOLUTION);
    options.min_spatial_resolution(MIN_SPATIAL_RESOLUTION);
    options.temporal_resolution(TEMPORAL_RESOLUTION);
    const Indexer::Options current = sbs.indexer_options();

    Logger::info() << "Requested database configuration update:"
                   << "\n  Maximum spatial resolution (deg): " << options.max_spatial_resolution() / DEG
                   << " -> " << current.max_spatial_resolution() / DEG
                   << "\n  Minimum spatial resolution (deg): " << options.min_spatial_resolution() / DEG
                   << " -> " << current.min_spatial_resolution() / DEG
                   << "\n  Temporal resolution (1/day): " << options.temporal_resolution()
                   << " -> " << current.temporal_resolution()
                   << std::endl;

    if (options == current)
        Logger::info() << "Database configuration matches, no need to re-index." << std::endl;
    else
        sbs.reindex(options);

    return 0;
}