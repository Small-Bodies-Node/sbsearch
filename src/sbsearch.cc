#include "config.h"

#include <iostream>
#include <unordered_set>
#include <string>
#include <vector>
#include <s2/s2latlng.h>
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
            db = new SBSearchDatabaseSqlite3(name);

        db->setup_tables();
        indexer_ = Indexer(indexer_options);
    }

    void SBSearch::reindex()
    {
        int n;
        int64 i = 0;
        vector<int64> observation_ids;

        n = db->get_int64("SELECT COUNT(*) FROM observations");
        while (i < n)
        {
            observation_ids.clear();
            for (int j = 0; j < 1000; j++)
            {
                if (i == n)
                    break;

                observation_ids.push_back(++i);
            }
            vector<Observation> observations = db->get_observations(observation_ids.begin(), observation_ids.end());
        }
    }

    void SBSearch::add_observations(vector<Observation> &observations)
    {
        // index observations, as needed
        for (int i = 0; i < observations.size(); i++)
            if (observations[i].terms().size() == 0)
                observations[i].terms(indexer_.index_terms(observations[i].as_polygon(), observations[i].mjd_start(), observations[i].mjd_stop()));

        db->add_observations(observations);
    }

    std::pair<double, double> SBSearch::date_range(string source)
    {
        return db->date_range(source);
    }

    // Only searches the database by spatial index (not spatial-temporal).
    vector<Observation> SBSearch::find_observations(const S2Point &point, double mjd_start, double mjd_stop)
    {
        if ((mjd_start > mjd_stop) && (mjd_stop != -1))
            throw std::runtime_error("Temporal search requested, but mjd_start > mjd_stop.");

        vector<string> query_terms = indexer_.query_terms(point);
        vector<Observation> approximate_matches = db->find_observations(query_terms);

        // collect observations that cover point and are within the requested time range
        vector<Observation> matches;
        for (auto observation : approximate_matches)
        {
            // check dates, if requested
            if ((mjd_start != -1) & (observation.mjd_stop() < mjd_start))
                continue;

            if ((mjd_stop != -1) & ((observation.mjd_start() > mjd_stop)))
                continue;

            // check spatial intersection
            if (observation.as_polygon().Contains(point))
                matches.push_back(observation);
        }

        std::cout << "  Matched " << matches.size() << " of " << approximate_matches.size() << " approximate matches." << std::endl;

        return matches;
    }

    // Only searches the database by spatial index (not spatial-temporal).
    vector<Observation> SBSearch::find_observations(const S2Polygon &polygon, double mjd_start, double mjd_stop)
    {
        if (mjd_start > mjd_stop)
            throw std::runtime_error("Temporal search requested, but mjd_start > mjd_stop.");

        vector<string> query_terms = indexer_.query_terms(polygon);
        vector<Observation> approximate_matches = db->find_observations(query_terms);

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

    // Searches the database by spatial-temporal index.
    vector<Found> SBSearch::find_observations(const Ephemeris &eph)
    {
        std::vector<Observation> all_matches;
        for (auto segment : eph.segments())
        {
            vector<string> query_terms = indexer_.query_terms(segment);
            vector<Observation> segment_matches = db->find_observations(query_terms);
            all_matches.insert(all_matches.end(), segment_matches.begin(), segment_matches.end());
        }

        // Because the search is segment by segment, duplicates can accumulate.
        std::unordered_set<Observation> unique_matches(all_matches.begin(), all_matches.end());

        vector<Found> found;
        for (auto observation : unique_matches)
        {
            Ephemeris e = eph.subsample(observation.mjd_start(), observation.mjd_stop());
            if (observation.as_polygon().Intersects(e.as_polygon()))
                found.emplace_back(observation, e);
        }

        std::cout << "  Matched " << found.size() << " of " << unique_matches.size() << " approximate matches." << std::endl;
        return found;
    }
}