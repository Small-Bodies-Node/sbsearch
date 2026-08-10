#ifndef SBSEARCH_UTIL_STRING_H_
#define SBSEARCH_UTIL_STRING_H_

#include <charconv>
#include <string>
#include <string_view>
#include <vector>
#include <s2/s2latlng.h>
#include <s2/s2latlng_rect.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>

using std::string;
using std::string_view;
using std::vector;

namespace sbsearch::util
{
    // Split a string given delimiter.  The delimiter is not included in the
    // output.
    vector<string> split(std::string_view str, const char delimiter);

    // Strip leading and trailing spaces.
    string strip(std::string_view str);

    // Convert string_view to double;
    double svtod(string_view s);

    // Convert string_view to long;
    // double svtod(string_view s);

    // Join a vector with the delimiter.
    template <typename T>
    string join(const vector<T> &v, string_view delimiter);

    // String-formatted vertices, comma-separated RA:Dec pairs in units of
    // degrees, e.g., "0:0, 0:1, 1:1".  For a polygon, only the first loop is
    // used.
    string format_vertices(const vector<S2LatLng> &vertices);
    string format_vertices(const vector<S2Point> &vertices);
    string format_vertices(const S2LatLngRect &fov);
    string format_vertices(const S2Polygon &polygon);

    // Convert string format ("RA:Dec, ...", units of degrees) to vector of points
    vector<S2Point> make_vertices(string_view str);

    ////////////////////////////////////////////////////////////////////////////////
    // implementations
    template <typename T>
    string join(const vector<T> &v, string_view delimiter)
    {
        if (v.empty())
            return "";

        std::stringstream s;
        s << v.front();
        for (auto it = std::next(v.begin()); it < v.end(); it = std::next(it))
            s << delimiter << *it;
        return s.str();
    }
}

#endif // SBSEARCH_UTIL_STRING_H_