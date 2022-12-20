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
    void SBSearchDatabase::drop_time_indices()
    {
        execute_sql("DROP INDEX IF EXISTS idx_observations_mjd_start;"
                    "DROP INDEX IF EXISTS idx_observations_mjd_stop;");
    }

    void SBSearchDatabase::add_observations(vector<Observation> &observations)
    {
        execute_sql("BEGIN TRANSACTION;");
        for (auto observation : observations)
            add_observation(observation);
        execute_sql("END TRANSACTION;");
    }
}