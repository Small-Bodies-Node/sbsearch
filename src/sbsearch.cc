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

namespace sbsearch
{
    struct Found
    {
        Observation observation;
        Ephemeris ephemeris;
    };

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
            std::cerr << "passed temporal, checking spatial\n";
            if (observation.as_polygon().Intersects(polygon))
                matches.push_back(observation);
        }

        std::cerr << "  Matched " << matches.size() << " of " << approximate_matches.size() << " approximate matches." << std::endl;

        return matches;
    }
}