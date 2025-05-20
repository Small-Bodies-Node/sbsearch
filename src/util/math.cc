#include "config.h"

#include <algorithm>
#include <optional>
#include <vector>

#include "math.h"

#include <iostream>
#include <iterator>

using std::optional;
using std::vector;

namespace sbsearch::util
{
    bool is_increasing(const vector<double> &v)
    {
        auto geq = std::greater_equal<double>();
        auto i = std::adjacent_find(v.begin(), v.end(), geq);
        return (i == v.end());
    }

    bool is_increasing(const vector<optional<double>> &vo)
    {
        try
        {
            vector<double> v;
            v.reserve(vo.size());
            std::transform(vo.begin(), vo.end(), std::back_inserter(v),
                           [](const auto &value)
                           { return value.value(); });
            return is_increasing(v);
        }
        catch (const std::bad_optional_access &e)
        {
        }
        return false;
    }

    std::optional<double> interp(const std::optional<double> a, const std::optional<double> b, const double frac)
    {
        if (!a || !b)
            return std::nullopt;

        return interp(a.value(), b.value(), frac);
    }
}