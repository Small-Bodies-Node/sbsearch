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

    void SBSearchDatabase::Options::max_spatial_cells(int n) { max_spatial_cells_ = n; };
    int SBSearchDatabase::Options::max_spatial_cells() const { return max_spatial_cells_; };
    int SBSearchDatabase::Options::max_spatial_level() const { return max_spatial_level_; };
    int SBSearchDatabase::Options::min_spatial_level() const { return min_spatial_level_; };

    void SBSearchDatabase::Options::max_spatial_resolution(double arcmin)
    {
        min_spatial_level_ = S2::kAvgAngleSpan.GetLevelForMaxValue(arcmin * ARCMIN);
    };
    double SBSearchDatabase::Options::max_spatial_resolution() const { return S2::kAvgAngleSpan.GetValue(min_spatial_level_) / ARCMIN; };

    void SBSearchDatabase::Options::min_spatial_resolution(double arcmin)
    {
        max_spatial_level_ = S2::kAvgAngleSpan.GetLevelForMaxValue(arcmin * ARCMIN);
    };
    double SBSearchDatabase::Options::min_spatial_resolution() const { return S2::kAvgAngleSpan.GetValue(max_spatial_level_) / ARCMIN; };

    // temporal resolution (nearest integer fraction between 0.01 and 1)
    void SBSearchDatabase::Options::temporal_resolution(double days)
    {
        time_terms_per_day_ = std::max(1, std::min(100, (int)std::round(1.0 / days)));
    };
    double SBSearchDatabase::Options::temporal_resolution() const { return 1.0 / time_terms_per_day_; };

    SBSearchDatabase::SBSearchDatabase(const Options &options) : options_(options)
    {
        S2RegionTermIndexer::Options region_term_indexer_options;
        region_term_indexer_options.set_min_level(options_.min_spatial_level());
        region_term_indexer_options.set_max_level(options_.max_spatial_level());
        region_term_indexer_options.set_max_cells(options_.max_spatial_cells());
        indexer = S2RegionTermIndexer(region_term_indexer_options);

        cout << "\nSpatial index limits " << options_.min_spatial_resolution() << " - " << options_.max_spatial_resolution()
             << " arcmin (levels " << options_.min_spatial_level() << " - " << options_.max_spatial_level()
             << "); time index resolution: " << 86400 * options_.temporal_resolution() << " s.\n\n";
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