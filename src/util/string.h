#ifndef SBSEARCH_UTIL_STRING_H_
#define SBSEARCH_UTIL_STRING_H_

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
    // Make a plural string.
    string pluralize(int count, string_view singular, string_view plural);

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