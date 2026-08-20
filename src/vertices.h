#ifndef SBS_VERTICES_H_
#define SBS_VERTICES_H_

#include <string>
#include <string_view>
#include <vector>
#include <s2/s2latlng.h>
#include <s2/s2latlng_rect.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>

#include "ephemeris/ephemeris.h"

using sbsearch::ephemeris::Ephemeris;
using std::string;
using std::string_view;
using std::vector;

namespace sbsearch
{
    // String-formatted vertices, comma-separated RA:Dec pairs in units of
    // degrees, e.g., "0:0, 0:1, 1:1".  For a polygon, only the first loop is
    // used.
    string format_vertices(const vector<S2LatLng> &vertices);
    string format_vertices(const vector<S2Point> &vertices);
    string format_vertices(const S2LatLngRect &fov);
    string format_vertices(const S2Polygon &polygon);

    // Convert string format ("RA:Dec, ...", units of degrees) to vector of points
    vector<S2Point> make_vertices(string_view str);

    // Convert Ephemeris::Datum to S2LatLng.
    inline S2LatLng make_latlng(const Ephemeris::Datum &d)
    {
        return S2LatLng::FromDegrees(d.dec, d.ra).Normalized();
    }

    // Convert Ephemeris::Datum to S2Point
    inline S2Point make_point(const Ephemeris::Datum &d)
    {
        return make_latlng(d).ToPoint();
    }
}

#endif