#include "config.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
#include <s2/s1chord_angle.h>
#include <s2/s2builderutil_snap_functions.h>
#include <s2/s2cap.h>
#include <s2/s2latlng.h>
#include <s2/s2lax_polygon_shape.h>
#include <s2/s2polygon.h>
#include <s2/s2point.h>
#include <s2/s2region.h>
#include <s2/mutable_s2shape_index.h>
#include <sys/stat.h>

#include "cli.h"
#include "ephemeris.h"
#include "exceptions.h"
#include "indexer.h"
#include "logging.h"
#include "observation.h"
#include "query_info.h"
#include "sbsearch.h"
#include "sbsdb/sbsdb.h"
#include "util/polygon.h"
#include "util/string.h"

using sbsearch::sbsdb::Postgresql;
using std::endl;

namespace sbsearch
{
    template <typename SBSDB>
    void SBSearch<SBSDB>::reindex(const Indexer::Options &options)
    {
        auto n = sbsdb::count::observations(&db_, 0, 10000);
        Logger::info() << "Re-indexing " << n << " observations." << endl;

        sbsdb::update::indexer_options(&db_, options);
        Logger::warning() << "Database configuration has been updated." << endl;
        indexer_ = Indexer(options);

        db_.drop_observations_indices();
        ProgressPercent widget(n);

        vector<int64_t> observation_ids;
        vector<vector<string>> observation_terms;

        const int chunk = 10000;
        observation_ids.reserve(chunk);
        observation_terms.reserve(chunk);

        int64_t offset = 0;
        do
        {
            observation_ids.clear();
            observation_terms.clear();

            S2Polygon polygon;
            // get the next chunk of ids and terms
            for (auto [observation_id, fov] : sbsdb::get::all_observations_fov(&db_, chunk, offset))
            {
                observation_ids.emplace_back(observation_id);
                util::make_polygon(util::make_vertices(fov), polygon);
                observation_terms.emplace_back(indexer_.terms(Indexer::index, polygon));
            }

            // update database terms
            sbsdb::update::observations(&db_, observation_ids, observation_terms);
            offset += observation_terms.size();
        } while (observation_terms.size() > 0);

        db_.create_observations_indices();
    }

    template <typename SBSDB>
    void SBSearch<SBSDB>::add_ephemeris(Ephemeris &eph)
    {
        if (sbsdb::count::ephemeris(&db_, eph.target(), eph.data(0).mjd.value(), eph.data(-1).mjd.value()) != 0)
            throw EphemerisError("data already present in database for target and date range: " +
                                 eph.target().to_string() + ", " +
                                 std::to_string(eph.data(0).mjd.value()) + ", " +
                                 std::to_string(eph.data(-1).mjd.value()));

        sbsdb::add::ephemeris(&db_, eph);
    }

    template <typename SBSDB>
    void SBSearch<SBSDB>::index_observations(Observations &observations)
    {
        for (auto observation = observations.begin(); observation < observations.end(); observation++)
        {
            if (observation->terms().size() == 0)
                observation->terms(indexer_.terms(Indexer::index, *observation));

            if (!observation->center())
            {
                vector<S2Point> points = util::make_vertices(observation->fov());
                const S2Point point = (points[0] / 4 + points[1] / 4 + points[2] / 4 + points[3] / 4).Normalize();
                observation->center(center_indexer_.GetIndexTerms(point, "")[0]);
            }
        }
    }

    template <typename SBSDB>
    void SBSearch<SBSDB>::add_observations(Observations &observations)
    {
        index_observations(observations);
        sbsdb::add::observations(&db_, observations);
    }

    template <typename SBSDB>
    void SBSearch<SBSDB>::update_observations(Observations &observations)
    {
        index_observations(observations);
        sbsdb::update::observations(&db_, observations);
    }

    template <typename SBSDB>
    Observations SBSearch<SBSDB>::find_observations(const S2Point &point, const FindOptions &options)
    {
        options.validate();

        if (options.padding > 0)
        {
            S2Cap cap(point, S1ChordAngle::Degrees(options.padding / 60));
            return find_observations(cap, options);
        }

        indexer_.mutable_options().max_spatial_query_cells(options.max_spatial_query_cells);

        vector<string> query_terms = indexer_.terms(Indexer::query, point);

        sbsdb::find::observations(&db_, query_terms, options.as_sbsearch_db_options());
        Observations matches = sbsdb::find::results(&db_);

        // only need approximate results?  done!
        if (options.approximate)
            return matches;

        int n_approximate_matches = matches.size();

        // keep observations that actually cover the point and are within the
        // requested time range
        auto not_intersecting = [&](const Observation &obs)
        { return !contains(obs, point, options.mjd_start, options.mjd_stop); };
        matches.data.erase(std::remove_if(matches.data.begin(), matches.data.end(), not_intersecting),
                           matches.data.end());

        Logger::info() << "Matched " << matches.size() << " of "
                       << n_approximate_matches << " approximate matches." << endl;

        return matches;
    }

    template <typename SBSDB>
    Observations SBSearch<SBSDB>::find_observations(const S2Cap &cap, const FindOptions &options)
    {
        options.validate();

        indexer_.mutable_options().max_spatial_query_cells(options.max_spatial_query_cells);

        vector<string> query_terms = indexer_.terms(Indexer::query, cap);

        sbsdb::find::observations(&db_, query_terms, options.as_sbsearch_db_options());
        Observations matches = sbsdb::find::results(&db_);

        // only need approximate results?  done!
        if (options.approximate)
            return matches;

        int n_approximate_matches = matches.size();

        // keep observations that intersect the area and are within the
        // requested time range
        auto not_intersecting = [&](const Observation &obs)
        { return !intersects(obs, cap, options.intersection_type, options.mjd_start, options.mjd_stop); };
        matches.data.erase(std::remove_if(matches.data.begin(), matches.data.end(), not_intersecting),
                           matches.data.end());

        Logger::info() << "Matched " << matches.size() << " of "
                       << n_approximate_matches << " approximate matches." << endl;

        return matches;
    }

    template <typename SBSDB>
    Observations SBSearch<SBSDB>::find_observations(const S2Polygon &polygon, const FindOptions &options)
    {
        options.validate();

        indexer_.mutable_options().max_spatial_query_cells(options.max_spatial_query_cells);

        S2Polygon query_polygon;
        util::padded_polygon(polygon, options.padding, query_polygon);

        vector<string> query_terms = indexer_.terms(Indexer::query, query_polygon);
        sbsdb::find::observations(&db_, query_terms, options.as_sbsearch_db_options());
        Observations matches = sbsdb::find::results(&db_);

        // only need approximate results?  done!
        if (options.approximate)
            return matches;

        int n_approximate_matches = matches.size();

        // keep observations that intersect the area and are within the
        // requested time range
        auto not_intersecting = [&](const Observation &obs)
        { return !intersects(obs, query_polygon, options.intersection_type, options.mjd_start, options.mjd_stop); };

        matches.data.erase(std::remove_if(matches.data.begin(), matches.data.end(), not_intersecting),
                           matches.data.end());

        cli::message("Matched " + std::to_string(matches.size()) + " of " +
                     std::to_string(n_approximate_matches) + " approximate matches.");

        return matches;
    }

    template <typename SBSDB>
    Founds SBSearch<SBSDB>::find_observations(const Ephemeris &ephemeris, const FindOptions &options)
    {
        options.validate();

        // reset query info
        query_info_.reset();

        Observatories observatories = sbsdb::get::all_observatories(&db_);

        cli::message(
            "Searching for observations with ephemeris: " +
            std::to_string(ephemeris.as_polyline().GetLength().degrees()) + " deg, " +
            std::to_string(ephemeris.data(-1).mjd.value() - ephemeris.data(0).mjd.value()) + " days.");

        indexer_.mutable_options().max_spatial_query_cells(options.max_spatial_query_cells);

        // split ephemeris into search segments
        vector<Ephemeris> segments = ephemeris.split(options.arc_length, options.time_period);
        string message = "Ephemeris split into " + std::to_string(segments.size()) + " segments.";
        std::cout << message << endl;
        Logger::debug() << message << endl;

        // search for each segment
        std::set<string> query_terms;
        ProgressPercent progress(segments.size());
        for (auto const &segment : segments)
        {
            // Account for padding and possibly parallax.
            double padding = options.padding;
            if (options.parallax)
            {
                // Increase search area by the size of the Earth at the distance
                // of the target = 8.7" / Delta = 0.145' / Delta, for Delta in au.
                auto delta = segment.delta();
                padding += 0.145 / std::max_element(delta.begin(), delta.end())->value();
            }

            vector<string> segment_query_terms = indexer_.terms(Indexer::query, segment, padding);

            auto db_options = options.as_sbsearch_db_options();
            db_options.mjd_start = segment.data(0).mjd.value();
            db_options.mjd_stop = segment.data(-1).mjd.value();
            sbsdb::find::observations(&db_, segment_query_terms, db_options);

            save_polygons(std::move(segment.as_polygons(padding)), options);
            save_terms(segment_query_terms, query_info_.query_terms, options);
            save_ephemeris(segment, options);

            progress.update();
            progress.status(false);
        }
        std::cout << "\n";

        Observations matches = sbsdb::find::results(&db_);
        save_terms(matches, query_info_.approximate_matches_index_terms, options);
        Logger::debug() << matches.size() << " approximate matches." << endl;

        // check for detailed intersection between ephemeris and candidates,
        // accounting for parallax in detail
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

            // only need approximate results?  done!
            if (options.approximate)
            {
                founds.append(Found(observation, eph));
                continue;
            }

            if (intersects(observation, eph, options.padding, options.mjd_start, options.mjd_stop))
                founds.data.emplace_back(observation, eph);
        }
        save_terms(founds, query_info_.matches_index_terms, options);

        cli::message("Matched " + std::to_string(founds.size()) + " of " +
                     std::to_string(matches.size()) + " approximate matches.");

        if (options.save)
        {
            sbsdb::add::found(&db_, founds);
            Logger::info() << matches.size() << " found observations saved to the database." << endl;
        }

        return founds;
    }

    template <typename SBSDB>
    const QueryInfo SBSearch<SBSDB>::query_info()
    {
        return query_info_;
    }

    template <typename SBSDB>
    void SBSearch<SBSDB>::save_polygon(const std::unique_ptr<S2Polygon> &polygon, const FindOptions &options)
    {
        if (!options.save_info)
            return;

        array<array<double, 2>, 4> vertices;
        for (int i = 0; i < 4; i++)
            vertices[i] = {S2LatLng::Longitude(polygon->loop(0)->vertex(i)).degrees(),
                           S2LatLng::Latitude(polygon->loop(0)->vertex(i)).degrees()};
        query_info_.query_polygons.emplace_back(vertices);
    }

    template <typename SBSDB>
    void SBSearch<SBSDB>::save_polygons(const vector<std::unique_ptr<S2Polygon>> &polygons, const FindOptions &options)
    {
        if (!options.save_info)
            return;

        for (auto &polygon : polygons)
        {
            save_polygon(polygon, options);
        }
    }

    template <typename SBSDB>
    void SBSearch<SBSDB>::save_terms(const vector<string> &terms, set<string> &dest, const FindOptions &options)
    {
        if (!options.save_info)
            return;

        std::copy(terms.begin(), terms.end(), std::inserter(dest, dest.end()));
    }

    template <typename SBSDB>
    void SBSearch<SBSDB>::save_terms(const Observations &observations, set<string> &dest, const FindOptions &options)
    {
        if (!options.save_info)
            return;

        for (auto const &observation : observations)
            save_terms(observation.terms(), dest, options);
    }

    template <typename SBSDB>
    void SBSearch<SBSDB>::save_terms(const Founds &founds, set<string> &dest, const FindOptions &options)
    {
        if (!options.save_info)
            return;

        for (auto const &found : founds)
            save_terms(found.observation.terms(), dest, options);
    }

    template <typename SBSDB>
    void SBSearch<SBSDB>::save_ephemeris(const Ephemeris &ephemeris, const FindOptions &options)
    {
        if (!options.save_info)
            return;

        vector<array<double, 5>> vector;
        for (auto const &eph : ephemeris.data())
            vector.emplace_back(
                array<double, 5>{eph.ra.value_or(1e99),
                                 eph.dec.value_or(1e99),
                                 eph.unc_a.value_or(1e99),
                                 eph.unc_b.value_or(1e99),
                                 eph.unc_theta.value_or(1e99)});

        query_info_.ephemeris_segments.emplace_back(std::move(vector));
    }

    template void SBSearch<Postgresql>::reindex(const Indexer::Options &);
    template void SBSearch<Postgresql>::add_ephemeris(Ephemeris &);
    template void SBSearch<Postgresql>::index_observations(Observations &);
    template void SBSearch<Postgresql>::add_observations(Observations &);
    template void SBSearch<Postgresql>::update_observations(Observations &);
    template Observations SBSearch<Postgresql>::find_observations(const S2Point &, const FindOptions &);
    template Observations SBSearch<Postgresql>::find_observations(const S2Cap &, const FindOptions &);
    template Observations SBSearch<Postgresql>::find_observations(const S2Polygon &, const FindOptions &);
    template Founds SBSearch<Postgresql>::find_observations(const Ephemeris &, const FindOptions &);
    template const QueryInfo SBSearch<Postgresql>::query_info();
    template void SBSearch<Postgresql>::save_polygon(const std::unique_ptr<S2Polygon> &, const FindOptions &);
    template void SBSearch<Postgresql>::save_polygons(const vector<std::unique_ptr<S2Polygon>> &, const FindOptions &);
    template void SBSearch<Postgresql>::save_terms(const vector<string> &, set<string> &, const FindOptions &);
    template void SBSearch<Postgresql>::save_terms(const Observations &, set<string> &, const FindOptions &);
    template void SBSearch<Postgresql>::save_terms(const Founds &, set<string> &, const FindOptions &);
    template void SBSearch<Postgresql>::save_ephemeris(const Ephemeris &, const FindOptions &);
}