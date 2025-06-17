#ifndef SBSEARCH_UTIL_STRING_H_
#define SBSEARCH_UTIL_STRING_H_

#include <vector>
#include <string>
#include <s2/s2latlng.h>
#include <s2/s2latlng_rect.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>

using std::string;
using std::vector;

namespace sbsearch::util
{
    // Split a string given delimiter.  The delimiter is not included in the
    // output.
    vector<string> split(std::string_view str, const char delimiter);

    // Strip leading and trailing spaces.
    string strip(std::string_view str);

    // Join a vector with the delimiter.
    template <typename T>
    string join(const vector<T> &v, const string &delimiter);

    // String-formatted vertices, comma-separated RA:Dec pairs in units of
    // degrees, e.g., "0:0, 0:1, 1:1".  For a polygon, only the first loop is
    // used.
    string format_vertices(const vector<S2LatLng> &vertices);
    string format_vertices(const vector<S2Point> &vertices);
    string format_vertices(const S2LatLngRect &fov);
    string format_vertices(const S2Polygon &polygon);

    // Convert string format ("RA:Dec, ...", units of degrees) to vector of points
    vector<S2Point> make_vertices(const string &str);

    // implementations
    template <typename T>
    string join(const vector<T> &v, const string &delimiter)
    {
        if (v.empty())
            return "";

        std::stringstream s;
        // if constexpr (std::is_same_v<T, string> == true)
        // {
        s << v.front();
        for (auto it = std::next(v.begin()); it < v.end(); it = std::next(it))
            s << delimiter << *it;
        // }
        // else
        // {
        //     s << std::to_string(v.front());
        //     for (auto it = std::next(v.begin()); it < v.end(); it = std::next(it))
        //         s += delimiter + std::to_string(*it);
        // }
        return s.str();
        // return std::accumulate(v.begin(), v.end(), std::begin(s),
        //                        [&delimiter](const string &a, const T &b)
        //                        { return a + delimiter + std::to_string(b); });

        // string s{""};
        // std::copy(v.begin(), v.end(), std::back_inserter(s),
        //           [&delimiter](const string &a, const T &b)
        //           { return a + delimiter + std::to_string(b); });
        // return s;
    }

    // Get and parse the next line of a CSV formatted stream.  Blank lines and
    // lines starting with # are skipped.  Leading and trailing whitespace is
    // stripped unless within a quoted string.  The maximum cell length is 1024
    // characters.  Returns the number of lines read and a vector of cells.
    vector<string> get_csv_cells(std::istream &is);

    // string join(const vector<string> &v, const string &delimiter)
    // {
    //     if (v.empty())
    //         return "";

    //     string s = v.front();
    //     for (auto it = std::next(v.begin()); it < v.end(); it = std::next(it))
    //         s += delimiter + *it;

    //     return s;
    // }

}

#endif // SBSEARCH_UTIL_STRING_H_