#include "config.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <s2/s2point.h>
#include <s2/s2cap.h>
#include <s2/s2metrics.h>
#include <s2/s2polyline.h>
#include <s2/s2polygon.h>
#include <s2/s2region_term_indexer.h>

#include "ephemeris.h"
#include "indexer.h"
#include "observation.h"
#include "sbsdb.h"
#include "util.h"

using std::cout;
using std::endl;
using std::string;
using std::unique_ptr;
using std::vector;

namespace sbsearch
{
    void SBSearchDatabase::drop_observations_indices()
    {
        execute_sql("DROP INDEX IF EXISTS idx_observations_mjd_start;"
                    "DROP INDEX IF EXISTS idx_observations_mjd_stop;");
    }

    void SBSearchDatabase::create_observations_indices()
    {
        // does not include the terms index
        execute_sql("CREATE INDEX IF NOT EXISTS idx_observations_mjd_start ON observations(mjd_start);\n"
                    "CREATE INDEX IF NOT EXISTS idx_observations_mjd_stop ON observations(mjd_stop);\n"
                    "CREATE UNIQUE INDEX IF NOT EXISTS idx_observations_source_product_id ON observations(source, product_id);\n"
                    "CREATE INDEX IF NOT EXISTS idx_observations_product_id ON observations(product_id);\n");
    }

    void SBSearchDatabase::add_observations(vector<Observation> &observations)
    {
        execute_sql("BEGIN TRANSACTION;");
        for (auto observation : observations)
            add_observation(observation);
        execute_sql("END TRANSACTION;");
    }

    Indexer::Options SBSearchDatabase::indexer_options()
    {
        Indexer::Options options;
        options.max_spatial_cells(get_int("SELECT value FROM configuration WHERE parameter=\"max_spatial_cells\";"));
        options.max_spatial_level(get_int("SELECT value FROM configuration WHERE parameter=\"max_spatial_level\";"));
        options.min_spatial_level(get_int("SELECT value FROM configuration WHERE parameter=\"min_spatial_level\";"));
        options.temporal_resolution(get_int("SELECT value FROM configuration WHERE parameter=\"temporal_resolution\";"));
        return options;
    }
}