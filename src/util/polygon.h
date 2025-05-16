#ifndef SBSEARCH_UTIL_POLYGON_H_
#define SBSEARCH_UTIL_POLYGON_H_

#include <vector>
#include <string>
#include <s2/s2point.h>
#include <s2/s2polygon.h>

using std::vector;

namespace sbsearch::util
{
    // Make a single simple polygon (no edge crossings, etc.).  The orientation
    // of the vertices will forced to make a polygon of area smaller than 2 pi.
    void make_polygon_simple(const vector<S2Point> &vertices, S2Polygon &polygon);

    // Make a polygon from a set of vertices will all edge and error checks.
    void make_polygon(const vector<S2Point> &vertices, S2Polygon &polygon);

    // Add a padding around the polygon given by pad in arcmin.  Padding must be
    // > 0 or else the original polygon will be returned unmodified.
    void padded_polygon(const S2Polygon &polygon, const double padding, S2Polygon &result);
}

#endif // SBSEARCH_UTIL_POLYGON_H_