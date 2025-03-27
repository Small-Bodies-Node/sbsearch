#include "config.h"

#include <algorithm>
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

#include "ephemeris.h"
#include "exceptions.h"
#include "logging.h"
#include "observation.h"
#include "sbsearch.h"
#include "sbsdb/sbsdb.h"
#include "util.h"

using std::endl;

namespace sbsearch
{
    std::istream &operator>>(std::istream &in, sbsearch::IntersectionType &intersection_type)
    {
        std::string token;
        in >> token;
        std::transform(token.begin(), token.end(), token.begin(),
                       [](unsigned char c)
                       { return std::tolower(c); });

        if ((token == "containspoint") || (token == "containscenter"))
            intersection_type = IntersectionType::ContainsPoint;
        else if (token == "containsarea")
            intersection_type = IntersectionType::ContainsArea;
        else if (token == "intersectsarea")
            intersection_type = IntersectionType::IntersectsArea;
        else if (token == "containedbyarea")
            intersection_type = IntersectionType::ContainedByArea;
        else
            in.setstate(std::ios_base::failbit);
        return in;
    }

    template <typename SBSDB>
    void SBSearch<SBSDB>::reindex(const Indexer::Options options)
    {
        auto n = sbsdb::count::observations(db_, 0, 10000);
        Logger::info() << "Re-indexing " << n << " observations." << endl;

        db_.drop_observations_indices();
        ProgressPercent widget(n);

        const int chunk = 10000;
        vector<vector<string>> terms;
        terms.reserve(chunk);

        int64_t offset = 0;
        do
        {
            terms.clear();

            S2Polygon polygon;
            for (auto [observation_id, fov] : sbsdb::get::all_observations_fov(db_, chunk, offset))
            {
                make_polygon(fov, polygon);
                terms.emplace_back(indexer_.terms(Indexer::index, polygon));
            }

            offset += terms.size();
        } while (terms.size() > 0);

        db_.create_observations_indices();

        // int64 i = 0;
        // const int chunk = 10000;
        // vector<int64_t> observation_ids;
        // observation_ids.reserve(chunk);

        // db_->indexer_options(options);
        // Logger::warning() << "Database configuration has been updated." << endl;
        // indexer_ = Indexer(options);

        // db_->drop_observations_indices();

        // ProgressPercent widget(n);
        // while (i < n)
        // {
        //     observation_ids.clear();
        //     for (int j = 0; j < chunk; j++)
        //     {
        //         if (i == n)
        //             break;

        //         observation_ids.push_back(++i);
        //     }
        //     Observations observations = db_->get_observations(observation_ids);

        //     // delete the terms and they will be regenerated
        //     for (vector<Observation>::iterator observation = observations.begin(); observation < observations.end(); observation++)
        //         observation->terms(vector<string>{});
        //     add_observations(observations);

        //     widget.update(observations.size());
        // }

        // db_->create_observations_indices();
    }

    template <typename SBSDB>
    void SBSearch<SBSDB>::add_ephemeris(Ephemeris &eph)
    {
        // MovingTarget target;

        // // validate ephemeris target
        // if (!eph.target().moving_target_id())
        // {
        //     // is the primary designation in the database?
        //     target = db_->get_moving_target(eph.target().designation());
        //     if (!target.moving_target_id())
        //     {
        //         // how about the alternate names?
        //         for (const string &name : target.alternate_names())
        //         {
        //             target = db_->get_moving_target(name);
        //             if (!target.moving_target_id())
        //             {
        //                 // found it, but warn the user that it required an alternate name
        //                 Logger::warning() << "Found target " << eph.target().designation()
        //                                   << " in the database under alternate name " << name
        //                                   << ".  Consider updating the database with a moving target merge operation."
        //                                   << std::endl;
        //                 break;
        //             }
        //         }
        //     }

        //     // if the target is still undefined, add it
        //     if (!target.moving_target_id())
        //     {
        //         eph.target(target);
        //         db_->add_moving_target(target);
        //     }
        // }
        // else
        // {
        //     // verify target is in the database by getting it with the moving target ID
        //     target = db_->get_moving_target(eph.target().moving_target_id().value());
        // }

        // update ephemeris object as needed
        // if (target != eph.target())
        //     eph.target(target);

        // verify that no ephemeris data is already in the database for this
        // date range and target
        // if (get::ephemeris(target, eph.data(0).mjd, eph.data(-1).mjd).num_vertices() != 0)

        if (sbsdb::count::ephemeris(db_, eph.target(), eph.data(0).mjd, eph.data(-1).mjd) != 0)
            throw EphemerisError("data already present in database for target and date range: " +
                                 eph.target().to_string() + ", " +
                                 std::to_string(eph.data(0).mjd) + ", " +
                                 std::to_string(eph.data(-1).mjd));

        sbsdb::add::ephemeris(db_, eph);
    }

    template <typename SBSDB>
    void SBSearch<SBSDB>::add_observations(Observations &observations)
    {
        // index observations, as needed
        for (vector<Observation>::iterator observation = observations.begin(); observation < observations.end(); observation++)
            if (observation->terms().size() == 0)
                observation->terms(indexer_.terms(Indexer::index, *observation));

        sbsdb::add::observations(db_, observations);
    }

    template <typename SBSDB>
    Observations SBSearch<SBSDB>::get_observations(const vector<int64_t> &observation_ids)
    {
        return sbsdb::get::observations(db_, observation_ids);
    }

    template <typename SBSDB>
    Observations SBSearch<SBSDB>::find_observations(const S2Point &point, const FindOptions &options)
    {
        S2Cap cap(point, S1ChordAngle::Degrees(options.padding / 60));

        // Only searches the database by spatial index (not spatial-temporal).
        if ((options.mjd_start > options.mjd_stop) && (options.mjd_stop != -1))
            throw std::runtime_error("Temporal search requested, but mjd_start > mjd_stop.");

        indexer_.mutable_options().max_spatial_query_cells(options.max_spatial_query_cells);

        vector<string> query_terms;
        if (options.padding <= 0)
            query_terms = indexer_.terms(Indexer::query, point);
        else
            query_terms = indexer_.terms(Indexer::query, cap);

        auto matches = sbsdb::find::observations(
            db_, query_terms, options.as_sbsearch_database_options());

        // keep observations that cover the point or intersect the area and are
        // within the requested time range
        auto not_intersecting = [&](auto const observation)
        {
            S2Polygon polygon;

            // check dates, if requested
            if ((options.mjd_start != -1) & (observation.mjd_stop() < options.mjd_start))
                return true;

            if ((options.mjd_stop != -1) & ((observation.mjd_start() > options.mjd_stop)))
                return true;

            // check spatial intersection
            observation.as_polygon(polygon);
            if (options.padding <= 0)
                return !polygon.Contains(point);
            else
                return !intersects(polygon, cap, options.intersection_type);
        };

        int approximate_matches = matches.size();
        const auto new_end = std::remove_if(matches.begin(), matches.end(), not_intersecting);
        matches.erase(new_end, matches.end());

        Logger::info() << "Matched " << matches.size() << " of "
                       << approximate_matches << " approximate matches." << endl;

        return {matches.begin(), matches.end()};
    }

    template <typename SBSDB>
    Observations SBSearch<SBSDB>::find_observations(const S2Polygon &polygon, const FindOptions &options)
    {
        // Only searches the database by spatial index (not spatial-temporal).
        if (options.mjd_start > options.mjd_stop)
            throw std::runtime_error("Temporal search requested, but mjd_start > mjd_stop.");

        indexer_.mutable_options().max_spatial_query_cells(options.max_spatial_query_cells);

        S2Polygon query_polygon;
        padded_polygon(polygon, options.padding, query_polygon);

        vector<string> query_terms = indexer_.terms(Indexer::query, query_polygon);
        auto matches = sbsdb::find::observations(
            db_, query_terms, options.as_sbsearch_database_options());

        // collect intersections
        S2Polygon fov_polygon;
        auto not_intersecting = [&](auto const &observation)
        {
            // check dates, if requested
            if ((options.mjd_start != -1) & (observation.mjd_stop() < options.mjd_start))
                return true;

            if ((options.mjd_stop != -1) & ((observation.mjd_start() > options.mjd_stop)))
                return true;

            // check detailed spatial intersection
            observation.as_polygon(fov_polygon);
            return !intersects(fov_polygon, query_polygon, options.intersection_type);
        };

        int approximate_matches = matches.size();
        const auto new_end = std::remove_if(matches.begin(), matches.end(), not_intersecting);
        matches.erase(new_end, matches.end());

        Logger::info() << "Matched " << matches.size() << " of "
                       << approximate_matches << " approximate matches." << endl;

        return {matches.begin(), matches.end()};
    }

    template <typename SBSDB>
    Founds SBSearch<SBSDB>::find_observations(const Ephemeris &ephemeris, const FindOptions &options)
    {
        Observatories observatories = sbsdb::get::all_observatories(db_);

        // Searches the database by spatial-temporal index.
        Logger::info() << "Searching for observations with ephemeris: "
                       << ephemeris.as_polyline().GetLength() << " deg, "
                       << (ephemeris.data(-1).mjd - ephemeris.data(0).mjd) << " days." << endl;

        indexer_.mutable_options().max_spatial_query_cells(options.max_spatial_query_cells);

        // split ephemeris into search segments
        vector<Ephemeris> segments = ephemeris.split(options.arc_length, options.time_period);
        Logger::debug() << "Ephemeris split into " << segments.size() << " segments." << endl;

        // search for each segment
        std::set<string> query_terms;
        std::unordered_set<Observation> matches;
        for (auto const &segment : segments)
        {
            // Account for padding and possibly parallax.
            double padding = options.padding;
            if (options.parallax)
            {
                // Increase search area by the size of the Earth at the distance
                // of the target = 8.7" / Delta = 0.145' / Delta, for Delta in au.
                auto delta = segment.delta();
                auto i = std::max_element(delta.begin(), delta.end());
                padding += 0.145 / *i;
            }

            vector<string> segment_query_terms = indexer_.terms(Indexer::query, segment, padding);

            auto db_options = options.as_sbsearch_database_options();
            db_options.mjd_start = segment.data(0).mjd();
            db_options.mjd_stop = segment.data(-1).mjd();
            matches.merge(sbsdb::find::observations(db_, segment_query_terms, db_options));
        }
        Logger::debug() << matches.size() << " approximate matches." << endl;

        // check for detailed intersection between ephemeris and candidates,
        // accounting for parallax in detail
        S2Polygon fov_polygon, segment_polygon, query_polygon;
        Founds founds;
        for (auto observation : matches)
        {
            // check dates
            if ((options.mjd_start != -1) & (observation.mjd_stop() < options.mjd_start))
                continue;

            if ((options.mjd_stop != -1) & ((observation.mjd_start() > options.mjd_stop)))
                continue;

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

            // If only approximate results are requested, then we are done.
            if (options.approximate)
            {
                founds.append(Found(observation, eph));
                continue;
            }

            // Otherwise, we check for detailed intersection.
            observation.as_polygon(fov_polygon);
            eph.as_polygon(segment_polygon);
            padded_polygon(segment_polygon, options.padding, query_polygon);

            if (intersects(fov_polygon, query_polygon, IntersectsArea))
                founds.append(Found(observation, eph));
        }

        Logger::info() << "Matched " << founds.size() << " of "
                       << matches.size() << " approximate matches." << endl;

        if (options.save)
        {
            sbsdb::add::found(db_, founds);
            Logger::info() << matches.size() << " found observations saved to the database." << endl;
        }

        return founds;
    }

    template <typename SBSDB>
    bool SBSearch<SBSDB>::intersects(const S2Polygon &polygon, const S2Cap &area, const IntersectionType intersection_type)
    {
        bool result = false;
        switch (intersection_type)
        {
        case (ContainsPoint):
            result = polygon.Contains(area.center());
            break;
        case (ContainsArea):
            result = (polygon.GetDistanceToBoundary(area.center()) > area.radius().ToAngle() & polygon.Contains(area.center()));
            break;
        case (IntersectsArea):
            result = polygon.GetDistance(area.center()) < area.radius().ToAngle();
            break;
        case (ContainedByArea):
            // only testing loop[0]; sbsearch does not use multiple loops
            const S2Loop *loop = polygon.loop(0);

            // check that each vertex is contained; immediately end loop if any
            // vertex is not
            result = true;
            for (int i = 0; i < loop->num_vertices(); i++)
            {
                result = area.InteriorContains(loop->vertex(i));
                if (!result)
                    break;
            }
            break;
        }
        return result;
    }

    template <typename SBSDB>
    bool SBSearch<SBSDB>::intersects(const S2Polygon &polygon, const S2Polygon &area, const IntersectionType intersection_type)
    {
        bool result = false;

        switch (intersection_type)
        {
        case (ContainsCenter):
            result = polygon.Contains(area.GetCentroid().Normalize());
            break;
        case (ContainsArea):
            result = polygon.Contains(area);
            break;
        case (IntersectsArea):
            result = polygon.Intersects(area);
            break;
        case (ContainedByArea):
            result = area.Contains(polygon);
            break;
        }
        return result;
    }
}