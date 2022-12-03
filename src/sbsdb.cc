#include "sbsdb.h"
#include "sbsearch.h"
#include "util.h"

#include "ephemeris.h"
#include "observation.h"

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

using std::cout;
using std::endl;
using std::string;
using std::unique_ptr;
using std::vector;

namespace sbsearch
{
    SBSearchDatabase::SBSearchDatabase()
    {
        S2RegionTermIndexer::Options options;
        options.set_max_level(S2::kAvgEdge.GetClosestLevel(SPATIAL_TERM_RESOLUTION * 0.00029089));
        options.set_max_cells(SPATIAL_INDEXER_MAX_CELLS);
        indexer = S2RegionTermIndexer(options);

        cout << "\nIndex setup:"
             << "\n  Spatial min level: " << options.min_level()
             << "\n  Spatial max level: " << options.max_level()
             << "\n  Time resolution: " << 24.0 / TIME_TERMS_PER_DAY << " hr\n\n";
    }

    void SBSearchDatabase::drop_time_indices()
    {
        execute_sql("DROP INDEX IF EXISTS idx_obs_mjdstart;"
                    "DROP INDEX IF EXISTS idx_obs_mjdstop;");
    }

    void SBSearchDatabase::add_observations(vector<Observation> &observations)
    {
        execute_sql("BEGIN TRANSACTION;");
        for (auto observation : observations)
            add_observation(observation);
        execute_sql("END TRANSACTION;");
    }

    vector<Found> SBSearchDatabase::find_observations(Ephemeris eph)
    {
        vector<Observation> observations = fuzzy_search(eph);
        vector<Found> found;
        S2Polyline polyline = eph.as_polyline();
        for (auto observation : observations)
        {
            unique_ptr<S2Polygon> polygon = observation.as_polygon();
            if (polygon->Intersects(polyline))
                found.push_back(Found{observation, eph});
        }
        return found;
    }

}