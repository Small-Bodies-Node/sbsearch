#include "config.h"

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
            if (!exists & !options.create)
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
        vector<int64> observation_ids;

        n = *(db_->get_int64("SELECT COUNT(*) FROM observations"));
        Logger::info() << "Re-indexing " << n << " observations." << endl;

        db_->indexer_options(options);
        Logger::warning() << "Database configuration has been updated." << endl;
        indexer_ = Indexer(options);

        ProgressPercent widget(n);
        while (i < n)
        {
            observation_ids.clear();
            for (int j = 0; j < 1000; j++)
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

        vector<string> query_terms = indexer_.query_terms(point);
        Observations approximate_matches = db_->find_observations(query_terms, options);

        // collect observations that cover point and are within the requested time range
        Observations matches;
        for (auto observation : approximate_matches)
        {
            // check dates, if requested
            if ((options.mjd_start != -1) & (observation.mjd_stop() < options.mjd_start))
                continue;

            if ((options.mjd_stop != -1) & ((observation.mjd_start() > options.mjd_stop)))
                continue;

            // check spatial intersection
            if (observation.as_polygon().Contains(point))
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

        vector<string> query_terms = indexer_.query_terms(polygon);
        Observations approximate_matches = db_->find_observations(query_terms, options);

        // collect intersections
        Observations matches;
        for (auto observation : approximate_matches)
        {
            // check dates
            if ((observation.mjd_start() > options.mjd_stop) | (observation.mjd_stop() < options.mjd_start))
                continue;

            // check detailed spatial intersection
            if (observation.as_polygon().Intersects(polygon))
                matches.push_back(observation);
        }

        Logger::info() << "Matched " << matches.size() << " of " << approximate_matches.size() << " approximate matches." << endl;

        return matches;
    }

    vector<Found> SBSearch::find_observations(const Ephemeris &ephemeris, const SearchOptions &options)
    {
        Observatories observatories = db_->get_observatories();

        // Searches the database by spatial-temporal index.
        Logger::info() << "Searching for observations with ephemeris: "
                       << ephemeris.as_polyline().GetLength() * DEG << " deg, "
                       << (ephemeris.data(-1).mjd - ephemeris.data(0).mjd) << " days." << std::endl;

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
        Observations matches = db_->find_observations(vector<string>(query_terms.begin(), query_terms.end()), options);

        vector<Found> found;
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

            if (observation.as_polygon().Intersects(eph.as_polygon()))
                found.emplace_back(observation, eph);
        }

        Logger::info() << "Matched " << found.size() << " of " << matches.size() << " approximate matches." << endl;
        return found;
    }

    std::ostream &operator<<(std::ostream &os, const Found &found)
    {
        // found.ephemeris is the segment that matches; interpolate it to observation mid-time
        double mjd = (found.observation.mjd_start() + found.observation.mjd_stop()) / 2;
        os << found.observation << "  " << found.ephemeris.interpolate(mjd);
        return os;
    }

    std::ostream &operator<<(std::ostream &os, const vector<Found> &founds)
    {
        // scan vector to determine column widths
        Observation::Format obs_format;
        int max_moving_target_id = 0;
        Ephemeris::Format eph_format;
        for (const Found &found : founds)
        {
            Observation::Format obsf = found.observation.format_widths();
            obs_format.source_width = std::max(obs_format.source_width, obsf.source_width);
            obs_format.observatory_width = std::max(obs_format.observatory_width, obsf.observatory_width);
            obs_format.product_id_width = std::max(obs_format.product_id_width, obsf.product_id_width);
            obs_format.fov_width = std::max(obs_format.fov_width, obsf.fov_width);
            obs_format.show_fov = std::max(obs_format.show_fov, obsf.show_fov);

            Ephemeris::Format ephf = found.ephemeris.format_widths();
            eph_format.designation_width = std::max(eph_format.designation_width, ephf.designation_width);
            eph_format.moving_target_id_width = std::max(eph_format.moving_target_id_width, ephf.moving_target_id_width);
            eph_format.tmtp_width = std::max(eph_format.tmtp_width, ephf.tmtp_width);
        }

        obs_format.source_width = std::max(obs_format.source_width, size_t(6));
        obs_format.observatory_width = std::max(obs_format.observatory_width, size_t(11));
        obs_format.observation_id_width = std::max(obs_format.observation_id_width, size_t(14));
        obs_format.product_id_width = std::max(obs_format.product_id_width, size_t(14));
        obs_format.exposure_time_width = std::max(obs_format.exposure_time_width, size_t(13));
        obs_format.quote_strings = false;

        eph_format.designation_width = (size_t)std::max(eph_format.designation_width, size_t(4));
        eph_format.moving_target_id_width = (size_t)std::max(eph_format.moving_target_id_width, size_t(16));

        // print headers
        os << std::setw(obs_format.observation_id_width)
           << "observation_id"
           << "  "
           << std::setw(obs_format.source_width)
           << "source"
           << "  "
           << std::setw(obs_format.observatory_width)
           << "observatory"
           << "  "
           << std::setw(obs_format.product_id_width)
           << "product_id"
           << "  "
           << std::setw(11)
           << "mjd_start"
           << "  "
           << std::setw(11)
           << "mjd_stop"
           << "  "
           << std::setw(obs_format.exposure_time_width)
           << "exposure_time"
           << "  ";

        if (obs_format.show_fov)
        {
            os << std::setw(obs_format.fov_width)
               << "fov"
               << "  ";
        }

        os << std::setw(eph_format.designation_width)
           << "desg"
           << "  "
           << std::setw(eph_format.moving_target_id_width)
           << "moving_target_id"
           << "  "
           << std::setw(11)
           << "mjd"
           << "  "
           << std::setw(eph_format.tmtp_width)
           << "tmtp"
           << "  "
           << std::setw(10)
           << "ra"
           << "  "
           << std::setw(10)
           << "dec"
           << "  "
           << std::setw(6)
           << "rh"
           << "  "
           << std::setw(6)
           << "delta"
           << "  "
           << std::setw(6)
           << "phase"
           << "\n"
           << std::setfill('-')
           << std::setw(obs_format.observation_id_width)
           << ""
           << "  "
           << std::setw(obs_format.source_width)
           << ""
           << "  "
           << std::setw(obs_format.observatory_width)
           << ""
           << "  "
           << std::setw(obs_format.product_id_width)
           << ""
           << "  "
           << std::setw(11)
           << ""
           << "  "
           << std::setw(11)
           << ""
           << "  "
           << std::setw(obs_format.exposure_time_width)
           << ""
           << "  ";

        if (obs_format.show_fov)
        {
            os << std::setw(obs_format.fov_width)
               << ""
               << "  ";
        }

        os << std::setw(eph_format.designation_width)
           << ""
           << "  "
           << std::setw(eph_format.moving_target_id_width)
           << ""
           << "  "
           << std::setw(11)
           << ""
           << "  "
           << std::setw(eph_format.tmtp_width)
           << ""
           << "  "
           << std::setw(10)
           << ""
           << "  "
           << std::setw(10)
           << ""
           << "  "
           << std::setw(6)
           << ""
           << "  "
           << std::setw(6)
           << ""
           << "  "
           << std::setw(6)
           << ""
           << "\n"
           << std::setfill(' ');

        for (Found found : founds)
        {
            found.observation.format = obs_format;
            found.ephemeris.format = eph_format;
            os << found << "\n";
        }
        return os;
    }
}