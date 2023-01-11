#include "config.h"

#include <iostream>
#include <string>
#include <vector>
#include <s2/s2polygon.h>
#include <s2/s2region.h>
#include <s2/mutable_s2shape_index.h>

#include "ephemeris.h"
#include "observation.h"
#include "sbsearch.h"
#include "sbsdb_sqlite3.h"

// debugging
#include "sbsearch_testing.h"

namespace sbsearch
{
    SBSearch::SBSearch(DatabaseType database_type, const char *name, Indexer::Options indexer_options)
    {
        if (database_type == sqlite3)
            db_ = new SBSearchDatabaseSqlite3(name);

        db_->setup_tables();
        indexer_ = Indexer(indexer_options);
    }

    void SBSearch::add_observations(vector<Observation> &observations)
    {
        // index observations, as needed
        for (int i = 0; i < observations.size(); i++)
            if (observations[i].terms().size() == 0)
                observations[i].terms(indexer_.index_terms(observations[i].as_polygon()));

        db_->add_observations(observations);
    }

    vector<Observation> SBSearch::find_observations(const S2Polygon &polygon, double mjd_start, double mjd_stop)
    {
        if (mjd_start > mjd_stop)
            throw std::runtime_error("Temporal search requested, but mjd_start > mjd_stop.");

        vector<string> query_terms = indexer_.query_terms(polygon);
        vector<Observation> approximate_matches = db_->find_observations(query_terms);

        // collect intersections
        vector<Observation> matches;
        for (auto observation : approximate_matches)
        {
            // check dates, if requested
            if ((mjd_stop != -1) & ((observation.mjd_start() > mjd_stop) | (observation.mjd_stop() < mjd_start)))
                continue;

            // check spatial intersection
            if (observation.as_polygon().Intersects(polygon))
                matches.push_back(observation);
        }

        std::cout << "  Matched " << matches.size() << " of " << approximate_matches.size() << " approximate matches." << std::endl;

        return matches;
    }

    vector<Found> SBSearch::find_observations(const Ephemeris &eph)
    {
        vector<Observation> approximate_matches;
        for (auto segment : eph.segments())
        {
            vector<string> query_terms = indexer_.query_terms(eph);
            vector<Observation> segment_approximate_matches = db_->find_observations(query_terms);
            approximate_matches.insert(approximate_matches.end(),
                                       segment_approximate_matches.begin(),
                                       segment_approximate_matches.end());

            std::cerr << segment.mjd(0) << " " << segment.mjd(1) << "\n";
            std::cerr << S2LatLng(segment.vertex(0)) << " " << S2LatLng(segment.vertex(1)) << "\n";
            testing::print_vector("terms ", query_terms);
        }

        vector<Found> found;
        for (auto observation : approximate_matches)
        {
            Ephemeris e = eph.subsample(observation.mjd_start(), observation.mjd_stop());
            if (observation.as_polygon().Intersects(e.as_polygon()))
                found.emplace_back(observation, e);

            std::cerr << e.mjd(0) << " " << e.mjd(1) << "\n";
            testing::print_vector("vertices ", e.vertices(), "\n");
        }

        std::cout << "  Matched " << found.size() << " of " << approximate_matches.size() << " approximate matches." << std::endl;
        return found;
    }
}