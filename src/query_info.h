#ifndef QUERYINFO_H_
#define QUERYINFO_H_

#include <array>
#include <set>
#include <string>
#include <vector>
#include <boost/json.hpp>

using std::array;
using std::set;
using std::string;
using std::vector;

namespace sbsearch
{
    // Information that may be used to debug or visualize an SBSearch query.
    struct QueryInfo
    {
        set<string> query_terms;
        set<string> approximate_matches_index_terms;
        set<string> matches_index_terms;

        // ra/dec, deg
        vector<array<array<double, 2>, 4>> query_polygons;

        // ra/dec/a/b/theta, deg/deg/arcsec/arcsec/deg
        vector<vector<array<double, 5>>> ephemeris_segments;

        void reset();
    };

    std::ostream &operator<<(std::ostream &os, const sbsearch::QueryInfo &info);
}

#endif