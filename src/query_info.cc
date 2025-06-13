#include "config.h"

#include <iostream>
#include <utility>
#include <vector>
#include <boost/json.hpp>

#include "query_info.h"

namespace sbsearch
{
    void QueryInfo::reset()
    {
        query_terms.erase(query_terms.begin(), query_terms.end());
        approximate_matches_index_terms.erase(approximate_matches_index_terms.begin(),
                                              approximate_matches_index_terms.end());
        matches_index_terms.erase(matches_index_terms.begin(),
                                  matches_index_terms.end());
        query_polygons.erase(query_polygons.begin(), query_polygons.end());
        ephemeris_segments.erase(ephemeris_segments.begin(), ephemeris_segments.end());
    }

    std::ostream &operator<<(std::ostream &os, const QueryInfo &info)
    {
        boost::json::object data;
        data["query terms"] = boost::json::array(
            info.query_terms.begin(),
            info.query_terms.end());
        data["approximate matches index terms"] = boost::json::array(
            info.approximate_matches_index_terms.begin(),
            info.approximate_matches_index_terms.end());
        data["matches index terms"] = boost::json::array(
            info.matches_index_terms.begin(),
            info.matches_index_terms.end());

        data["query polygons"].emplace_array();
        for (auto const &polygon : info.query_polygons)
        {
            boost::json::array vertices;
            for (auto const &[ra, dec] : polygon)
                vertices.push_back({ra, dec});

            data["query polygons"].as_array().push_back(vertices);
        }

        data["ephemeris segments"].emplace_array();
        for (auto const &segment : info.ephemeris_segments)
        {
            boost::json::array segment_data;
            for (auto const &[ra, dec, a, b, theta] : segment)
                segment_data.push_back({ra, dec, a, b, theta});
            data["ephemeris segments"].as_array().push_back(segment_data);
        }

        os << data;
        return os;
    }
}
