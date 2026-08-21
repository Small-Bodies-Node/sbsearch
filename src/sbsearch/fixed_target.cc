#include <algorithm>
#include <string>
#include <vector>
#include <set>
#include <s2/s1chord_angle.h>
#include <s2/s2cap.h>
#include <s2/s2latlng.h>
#include <s2/s2polygon.h>
#include <s2/s2point.h>

#include "cli.h"
#include "indexer.h"
#include "logging.h"
#include "observation.h"
#include "polygons.h"
#include "query_info.h"
#include "sbsdb.h"
#include "sbsearch.h"

using sbsearch::sbsdb::Postgresql;
using std::endl;

namespace sbsearch
{
    template <class SBSDB>
    Observations SBSearch<SBSDB>::search_observations(const S2Point &point, const SearchOptions &options)
    {
        options.validate();

        if (options.padding > 0)
        {
            S2Cap cap(point, S1ChordAngle::Degrees(options.padding / 60));
            return search_observations(cap, options);
        }

        indexer_.mutable_options().max_spatial_query_cells(options.max_spatial_query_cells);

        vector<string> query_terms = indexer_.terms(Indexer::query, point);

        sbsdb::search::observations(&db_, query_terms, options.as_sbsearch_db_options());
        Observations matches = sbsdb::search::results(&db_);

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

    template <class SBSDB>
    Observations SBSearch<SBSDB>::search_observations(const S2Cap &cap, const SearchOptions &options)
    {
        options.validate();

        indexer_.mutable_options().max_spatial_query_cells(options.max_spatial_query_cells);

        vector<string> query_terms = indexer_.terms(Indexer::query, cap);

        sbsdb::search::observations(&db_, query_terms, options.as_sbsearch_db_options());
        Observations matches = sbsdb::search::results(&db_);

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

    template <class SBSDB>
    Observations SBSearch<SBSDB>::search_observations(const S2Polygon &polygon, const SearchOptions &options)
    {
        options.validate();

        indexer_.mutable_options().max_spatial_query_cells(options.max_spatial_query_cells);

        S2Polygon query_polygon;
        padded_polygon(polygon, options.padding, query_polygon);

        vector<string> query_terms = indexer_.terms(Indexer::query, query_polygon);
        sbsdb::search::observations(&db_, query_terms, options.as_sbsearch_db_options());
        Observations matches = sbsdb::search::results(&db_);

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

        cli::message::debug("Matched " + std::to_string(matches.size()) + " of " +
                            std::to_string(n_approximate_matches) + " approximate matches.");

        return matches;
    }

    template Observations SBSearch<Postgresql>::search_observations(const S2Point &, const SearchOptions &);
    template Observations SBSearch<Postgresql>::search_observations(const S2Cap &, const SearchOptions &);
    template Observations SBSearch<Postgresql>::search_observations(const S2Polygon &, const SearchOptions &);
}