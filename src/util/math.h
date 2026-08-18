#ifndef SBSEARCH_UTIL_MATH_H_
#define SBSEARCH_UTIL_MATH_H_

#include <array>
#include <optional>
#include <vector>
#include <gsl/gsl_interp.h>
#include <s2/s2point.h>

using std::array;
using std::optional;
using std::vector;

namespace sbsearch::util
{
    // vector values must be increasing or constant
    bool is_monotonically_increasing(const vector<double> &v);

    // vector values must all be defined and increasing or constant
    bool is_monotonically_increasing(const vector<optional<double>> &v);

    // vector values must be increasing
    bool is_strictly_increasing(const vector<double> &v);

    // vector values must all be defined and increasing
    bool is_strictly_increasing(const vector<optional<double>> &v);

    // linear interpolation
    // there are no limits on `frac`, use `frac` < 0 or > 1 to extrapolate.
    inline double interp(const double a, const double b, const double frac)
    {
        return a + (b - a) * frac;
    };
    optional<double> interp(const optional<double> a, const optional<double> b, const double frac);

    // Check for and fix crossing edges given an arbitrary polygon with four
    // vertices.
    void fix_crossing_edges(vector<S2Point> &vertices);
}

#endif // SBSEARCH_UTIL_MATH_H_