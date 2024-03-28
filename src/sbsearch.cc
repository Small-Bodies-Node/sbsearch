#include "config.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
#include <s2/s2latlng.h>
#include <s2/s2polygon.h>
#include <s2/s2region.h>
#include <s2/mutable_s2shape_index.h>
#include <sys/stat.h>

#include "ephemeris.h"
#include "exceptions.h"
#include "logging.h"
#include "observation.h"
#include "sbsearch.h"
#include "sbsdb_sqlite3.h"

using std::endl;

namespace sbsearch
{
    SBSearch::SBSearch(DatabaseType database_type, const std::string name, const Options options)
    {
        // attempt to initialize logger
        Logger::get_logger(options.log_file).log_level(options.log_level);

        if (database_type == sqlite3)
        {
            struct stat buf;
            bool exists = (stat(name.c_str(), &buf) == 0);
            if (!exists & !options.create & name != ":memory:")
                throw std::runtime_error(name + " does not exist.");
            db_ = new SBSearchDatabaseSqlite3(name);
        }

        Indexer::Options indexer_options = db_->indexer_options();
        indexer_ = Indexer(indexer_options);
    }

    void SBSearch::reindex(const Indexer::Options options)
    {
        int n;
        int64 i = 0;
        const int N = 10000;
        vector<int64> observation_ids;
        observation_ids.reserve(N);

        n = *(db_->get_int64("SELECT COUNT(*) FROM observations"));
        Logger::info() << "Re-indexing " << n << " observations." << endl;

        db_->indexer_options(options);
        Logger::warning() << "Database configuration has been updated." << endl;
        indexer_ = Indexer(options);

        db_->drop_observations_indices();

        ProgressPercent widget(n);
        while (i < n)
        {
            observation_ids.clear();
            for (int j = 0; j < N; j++)
            {
                if (i == n)
                    break;

                observation_ids.push_back(++i);
            }
            Observations observations = db_->get_observations(observation_ids.begin(), observation_ids.end());

            // delete the terms and they will be regenerated
            for (Observation &observation : observations)
                observation.terms("");
            add_observations(observations);

            widget.update(observations.size());
        }

        db_->create_observations_indices();
    }

    void SBSearch::add_ephemeris(Ephemeris &eph)
    {
        MovingTarget target;

        // validate ephemeris target
        if (eph.target().moving_target_id() == UNDEF_MOVING_TARGET_ID)
        {
            // is the primary designation in the database?
            target = db_->get_moving_target(eph.target().designation());
            if (target.moving_target_id() == UNDEF_MOVING_TARGET_ID)
            {
                // how about the alternate names?
                for (const string &name : target.alternate_names())
                {
                    target = db_->get_moving_target(name);
                    if (target.moving_target_id() != UNDEF_MOVING_TARGET_ID)
                    {
                        // found it, but warn the user that it required an alternate name
                        Logger::warning() << "Found target " << eph.target().designation()
                                          << " in the database under alternate name " << name
                                          << ".  Consider updating the database with a moving target merge operation."
                                          << std::endl;
                        break;
                    }
                }
            }

            // if the target is still undefined, add it
            if (target.moving_target_id() == UNDEF_MOVING_TARGET_ID)
            {
                eph.target(target);
                db_->add_moving_target(target);
            }
        }
        else
        {
            // verify target is in the database by getting it with the moving target ID
            target = db_->get_moving_target(eph.target().moving_target_id());
        }

        // update ephemeris object as needed
        if (target != eph.target())
            eph.target(target);

        // verify that no ephemeris data is already in the database for this
        // date range and target
        if (db_->get_ephemeris(target, eph.data(0).mjd, eph.data(-1).mjd).num_vertices() != 0)
            throw EphemerisError("data already present in database for target and date range: " + target.designation() + ", " + std::to_string(eph.data(0).mjd) + ", " + std::to_string(eph.data(-1).mjd));

        db_->add_ephemeris(eph);
    }

    void SBSearch::add_observations(Observations &observations)
    {
        // index observations, as needed
        for (Observation &observation : observations)
            if (observation.terms().size() == 0)
                observation.terms(indexer_.index_terms(observation));

        db_->add_observations(observations);
    }

    Observations SBSearch::get_observations(const vector<int64> &observation_ids)
    {
        return db_->get_observations(observation_ids.begin(), observation_ids.end());
    }

    Observations SBSearch::find_observations(const S2Point &point, const SearchOptions &options)
    {
        // Only searches the database by spatial index (not spatial-temporal).
        if ((options.mjd_start > options.mjd_stop) && (options.mjd_stop != -1))
            throw std::runtime_error("Temporal search requested, but mjd_start > mjd_stop.");

        indexer_.mutable_options().max_spatial_query_cells(options.max_spatial_query_cells);

        vector<string> query_terms = indexer_.query_terms(point);
        Observations approximate_matches = db_->find_observations(query_terms, options.as_sbsearch_database_options());

        // collect observations that cover point and are within the requested time range
        S2Polygon polygon;
        Observations matches;
        for (auto observation : approximate_matches)
        {
            // check dates, if requested
            if ((options.mjd_start != -1) & (observation.mjd_stop() < options.mjd_start))
                continue;

            if ((options.mjd_stop != -1) & ((observation.mjd_start() > options.mjd_stop)))
                continue;

            // check spatial intersection
            observation.as_polygon(polygon);
            if (polygon.Contains(point))
                matches.push_back(observation);
        }

        Logger::info() << "Matched " << matches.size() << " of " << approximate_matches.size() << " approximate matches." << endl;

        return matches;
    }

    Observations SBSearch::find_observations(const S2Polygon &polygon, const SearchOptions &options)
    {
        // Only searches the database by spatial index (not spatial-temporal).
        if (options.mjd_start > options.mjd_stop)
            throw std::runtime_error("Temporal search requested, but mjd_start > mjd_stop.");

        indexer_.mutable_options().max_spatial_query_cells(options.max_spatial_query_cells);

        vector<string> query_terms = indexer_.query_terms(polygon);
        Observations approximate_matches = db_->find_observations(query_terms, options.as_sbsearch_database_options());

        // collect intersections
        S2Polygon fov_polygon;
        Observations matches;
        for (auto observation : approximate_matches)
        {
            // check dates
            if ((observation.mjd_start() > options.mjd_stop) | (observation.mjd_stop() < options.mjd_start))
                continue;

            // check detailed spatial intersection
            observation.as_polygon(fov_polygon);
            if (fov_polygon.Intersects(polygon))
                matches.push_back(observation);
        }

        Logger::info() << "Matched " << matches.size() << " of " << approximate_matches.size() << " approximate matches." << endl;

        return matches;
    }

    Founds SBSearch::find_observations(const Ephemeris &ephemeris, const SearchOptions &options)
    {
        Observatories observatories = db_->get_observatories();

        // Searches the database by spatial-temporal index.
        Logger::info() << "Searching for observations with ephemeris: "
                       << ephemeris.as_polyline().GetLength() << " deg, "
                       << (ephemeris.data(-1).mjd - ephemeris.data(0).mjd) << " days." << std::endl;

        indexer_.mutable_options().max_spatial_query_cells(options.max_spatial_query_cells);

        std::set<string> query_terms;
        for (auto segment : ephemeris.segments())
        {
            // Account for parallax?
            if (options.parallax)
            {
                // Increase search area by the size of the Earth at the distance
                // of the target = 8.7" / Delta, for Delta in au.
                const double delta_max = std::max({segment.data(0).delta, segment.data(1).delta});
                segment.mutable_options()->padding += 8.7 / delta_max;
            }
            vector<string> segment_query_terms = indexer_.query_terms(segment);
            query_terms.insert(segment_query_terms.begin(), segment_query_terms.end());
        }
        Observations matches = db_->find_observations(vector<string>(query_terms.begin(), query_terms.end()), options.as_sbsearch_database_options());

        S2Polygon fov_polygon, eph_polygon;
        Founds founds;
        for (auto observation : matches)
        {
            Ephemeris eph;

            // Approximate matches could have observation times beyond the
            // ephemeris bounds.  These observations should not be matched.  If
            // the user wanted such data, they should have requested a broader
            // ephemeris time span.
            try
            {
                eph = ephemeris.subsample(observation.mjd_start(), observation.mjd_stop());
            }
            catch (const std::runtime_error &)
            {
                continue;
            }

            // Account for parallax?  Then offset the ephemeris.
            if (options.parallax)
            {
                Observatory observatory;
                try
                {
                    observatory = observatories.at(observation.observatory());
                }
                catch (const std::out_of_range &)
                {
                    throw ObservatoryError(observation.observatory() + " not in database");
                }

                eph = eph.parallax_offset(observatory);
            }

            observation.as_polygon(fov_polygon);
            eph.as_polygon(eph_polygon);
            if (fov_polygon.Intersects(eph_polygon))
                founds.append(Found(observation, eph));
        }

        Logger::info() << "Matched " << founds.size() << " of " << matches.size() << " approximate matches." << endl;

        if (options.save)
        {
            db_->add_founds(founds);
            Logger::info() << matches.size() << " found observations saved to the database." << endl;
        }

        return founds;
    }
}