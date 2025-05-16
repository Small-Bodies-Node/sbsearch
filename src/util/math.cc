#include "config.h"

#include <algorithm>
#include <optional>
#include <vector>

#include "math.h"

using std::vector;

namespace sbsearch::util
{
    bool is_increasing(const vector<double> &v)
    {
        auto i = std::adjacent_find(v.begin(), v.end(), std::greater<double>());
        return (i == v.end());
    }

    double interp(const std::optional<double> a, const std::optional<double> b, const double frac)
    {
        if (!a || !b)
            return 0;

        return interp(a.value(), b.value(), frac);
    }
}