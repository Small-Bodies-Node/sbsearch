#include "config.h"

#include <iostream>

#include "indexer.h"
#include "sbsearch.h"
#include "test_db.h"

using sbsearch::Indexer;
using sbsearch::SBSearch;
using std::cout;

int main(int argc, char **argv)
{
    Indexer::Options options;
    options.max_spatial_cells(MAX_SPATIAL_CELLS);
    options.max_spatial_resolution(MAX_SPATIAL_RESOLUTION);
    options.min_spatial_resolution(MIN_SPATIAL_RESOLUTION);
    options.temporal_resolution(TEMPORAL_RESOLUTION);
    SBSearch sbs(SBSearch::sqlite3, "sbsearch_test.db", options);

    cout << "\nIndex setup:"
         << "\n  Minimum spatial resolution (deg): " << MIN_SPATIAL_RESOLUTION / DEG
         << "\n  Maximum spatial resolution (deg): " << MAX_SPATIAL_RESOLUTION / DEG
         << "\n  Temporal resolution (1/day): " << TEMPORAL_RESOLUTION
         << "\n\n";

    sbs.reindex();

    return 0;
}