#include "config.h"

#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>
#include <s2/s2latlng.h>
#include <s2/s2polygon.h>
#include <s2/s2region.h>
#include <s2/mutable_s2shape_index.h>

#include "ephemeris.h"
#include "logging.h"
#include "observation.h"
#include "sbsearch.h"
#include "sbsdb_sqlite3.h"

using std::endl;

namespace sbsearch
{
    SBSearch::SBSearch(DatabaseType database_type, const std::string name, const std::string log_file)
    {
        Logger::get_logger(log_file); // attempt to initialize logger

        if (database_type == sqlite3)
            db = new SBSearchDatabaseSqlite3(name);

        db->setup_tables();

        Indexer::Options options = db->indexer_options();
        indexer_ = Indexer(options);
    }

    void SBSearch::reindex(const Indexer::Options options)
    {
        int n;
        int64 i = 0;
        vector<int64> observation_ids;

        n = *(db->get_int64("SELECT COUNT(*) FROM observations"));
        Logger::info() << "Re-indexing " << n << " observations." << endl;

        db->indexer_options(options);
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
            vector<Observation> observations = db->get_observations(observation_ids.begin(), observation_ids.end());

            // delete the terms and they will be regenerated
            for (Observation &observation : observations)
                observation.terms("");
            add_observations(observations);

            widget.update(observations.size());
        }
    }

    void SBSearch::add_observations(vector<Observation> &observations)
    {
        // index observations, as needed
        for (int i = 0; i < observations.size(); i++)
            if (observations[i].terms().size() == 0)
            {
                observations[i].terms(indexer_.index_terms(observations[i].as_polygon(), observations[i].mjd_start(), observations[i].mjd_stop()));
            }

        db->add_observations(observations);
    }

    vector<Observation> SBSearch::get_observations(const vector<int64> &observation_ids)
    {
        return db->get_observations(observation_ids.begin(), observation_ids.end());
    }

    std::pair<double *, double *> SBSearch::date_range(string source)
    {
        return std::move(db->date_range(source));
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

        Logger::info() << "Matched " << matches.size() << " of " << approximate_matches.size() << " approximate matches." << endl;

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

        Logger::info() << "Matched " << matches.size() << " of " << approximate_matches.size() << " approximate matches." << endl;

        return matches;
    }

    // Searches the database by spatial-temporal index.
    vector<Found> SBSearch::find_observations(const Ephemeris &eph)
    {
        std::set<string> query_terms;
        for (auto segment : eph.segments())
        {
            vector<string> segment_query_terms = indexer_.query_terms(segment);
            query_terms.insert(segment_query_terms.begin(), segment_query_terms.end());
        }
        vector<Observation> matches = db->find_observations(vector<string>(query_terms.begin(), query_terms.end()));

        vector<Found> found;
        for (auto observation : matches)
        {
            Ephemeris e;

            // Approximate matches could have observation times beyond the
            // ephemeris bounds.  These observations should not be matched.  If
            // the user wanted such data, they should have requested a broader
            // ephemeris time span.
            try
            {
                e = eph.subsample(observation.mjd_start(), observation.mjd_stop());
            }
            catch (const std::runtime_error &)
            {
                continue;
            }

            if (observation.as_polygon().Intersects(e.as_polygon()))
                found.emplace_back(observation, e);
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
        int max_object_id = 0;
        for (const Found &found : founds)
        {
            Observation::Format _format = found.observation.format_widths();
            obs_format.source_width = std::max(obs_format.source_width, _format.source_width);
            obs_format.product_id_width = std::max(obs_format.product_id_width, _format.product_id_width);
            obs_format.fov_width = std::max(obs_format.fov_width, _format.fov_width);
            obs_format.show_fov = std::max(obs_format.show_fov, _format.show_fov);

            max_object_id = std::max(max_object_id, found.ephemeris.object_id());
        }
        obs_format.observation_id_width = std::max(obs_format.observation_id_width, size_t(14));
        obs_format.product_id_width = std::max(obs_format.product_id_width, size_t(14));
        obs_format.exposure_time_width = std::max(obs_format.exposure_time_width, size_t(13));
        obs_format.quote_strings = false;

        Ephemeris::Format eph_format = {
            std::max(size_t(std::floor(std::log10(max_object_id))) + 1, size_t(9))};

        // print headers
        os << std::setw(obs_format.observation_id_width)
           << "observation_id"
           << "  "
           << std::setw(obs_format.source_width)
           << "source"
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

        os << std::setw(eph_format.object_id_width)
           << "object_id"
           << "  "
           << std::setw(11)
           << "mjd"
           << "  "
           << std::setw(12)
           << "ra"
           << "  "
           << std::setw(12)
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

        os << std::setw(eph_format.object_id_width)
           << ""
           << "  "
           << std::setw(11)
           << ""
           << "  "
           << std::setw(12)
           << ""
           << "  "
           << std::setw(12)
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