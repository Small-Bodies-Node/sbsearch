#ifndef SBSEARCH_UTIL_MATH_H_
#define SBSEARCH_UTIL_MATH_H_

#include <vector>

using std::vector;

namespace sbsearch::util
{
    // vector values must always be increasing
    bool is_increasing(const vector<double> &v);

    // linear interpolation
    // there are no limits on `frac`, use `frac` < 0 or > 1 to extrapolate.
    inline double interp(const double a, const double b, const double frac)
    {
        return a + (b - a) * frac;
    };
    double interp(const std::optional<double> a, const std::optional<double> b, const double frac);
}

#endif // SBSEARCH_UTIL_MATH_H_