#include "config.h"

#include <algorithm>
#include <array>
#include <optional>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <s2/s2edge_crosser.h>

#include "util/math.h"
#include "util/optional.h"

#include <iostream>
#include <iterator>

using std::optional;
using std::vector;

namespace sbsearch::util
{
    bool is_monotonically_increasing(const vector<double> &v)
    {
        if (v.size() < 2)
            return true;

        auto test = std::greater<double>();
        auto i = std::adjacent_find(v.begin(), v.end(), test);
        return (i == v.end());
    }

    bool is_monotonically_increasing(const vector<optional<double>> &vo)
    {
        if (vo.size() < 2)
            return true;

        try
        {
            return is_monotonically_increasing(optionals_to_values(vo));
        }
        catch (const std::bad_optional_access &e)
        {
        }
        return false;
    }

    bool is_strictly_increasing(const vector<double> &v)
    {
        if (v.size() < 2)
            return true;

        auto test = std::greater_equal<double>();
        auto i = std::adjacent_find(v.begin(), v.end(), test);
        return (i == v.end());
    }

    bool is_strictly_increasing(const vector<optional<double>> &vo)
    {
        if (vo.size() < 2)
            return true;

        try
        {
            return is_strictly_increasing(optionals_to_values(vo));
        }
        catch (const std::bad_optional_access &e)
        {
        }
        return false;
    }

    optional<double> interp(const optional<double> a, const optional<double> b, const double frac)
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
            // swap 1/2
            std::iter_swap(vertices.begin() + 1, vertices.begin() + 2);
            return;
        }

        crosser.Init(&vertices[1], &vertices[2]);
        if (crosser.CrossingSign(&vertices[3], &vertices[0]) > 0)
        {
            // swap 2/3
            std::iter_swap(vertices.begin() + 2, vertices.begin() + 3);
            return;
        }
    }
}