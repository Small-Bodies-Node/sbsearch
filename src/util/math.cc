#include "config.h"

#include <algorithm>
#include <optional>
#include <vector>
#include <s2/s2edge_crosser.h>

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

    void fix_crossing_edges(vector<S2Point> &vertices)
    {
        // Check for crossing edges and re-order vertices as needed.
        S2EdgeCrosser crosser(&vertices[0], &vertices[1]);
        if (crosser.CrossingSign(&vertices[2], &vertices[3]) > 0)
        {
            std::cerr << "swapping 1/2\n";
            // swap 1/2
            std::iter_swap(vertices.begin() + 1, vertices.begin() + 2);
            return;
        }

        crosser.Init(&vertices[1], &vertices[2]);
        if (crosser.CrossingSign(&vertices[3], &vertices[0]) > 0)
        {
            std::cerr << "swapping 2/3\n";
            // swap 2/3
            std::iter_swap(vertices.begin() + 2, vertices.begin() + 3);
            return;
        }
    }
}