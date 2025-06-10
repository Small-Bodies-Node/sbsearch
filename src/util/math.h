#ifndef SBSEARCH_UTIL_MATH_H_
#define SBSEARCH_UTIL_MATH_H_

#include <optional>
#include <vector>
#include <s2/s2point.h>

using std::vector;

namespace sbsearch::util
{
    // vector values must always be increasing
    bool is_increasing(const vector<double> &v);

    // vector values must all be defined and always be increasing
    bool is_increasing(const vector<std::optional<double>> &v);

    // linear interpolation
    // there are no limits on `frac`, use `frac` < 0 or > 1 to extrapolate.
    inline double interp(const double a, const double b, const double frac)
    {
        return a + (b - a) * frac;
    };
    std::optional<double> interp(const std::optional<double> a, const std::optional<double> b, const double frac);

    // Check for and fix crossing edges given an arbitrary polygon with four
    // vertices.
    void fix_crossing_edges(vector<S2Point> &vertices);
}

#endif // SBSEARCH_UTIL_MATH_H_