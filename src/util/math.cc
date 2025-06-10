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

        // // Calculate the cross-product of the edges at each corner.
        // const S2Point A = vertices[1] - vertices[0];
        // const S2Point B = vertices[2] - vertices[1];
        // const S2Point C = vertices[3] - vertices[2];
        // const S2Point D = vertices[0] - vertices[3];
        // const S2Point corner0 = A.CrossProd(-D);
        // const S2Point corner1 = B.CrossProd(-A);
        // const S2Point corner2 = C.CrossProd(-B);
        // const S2Point corner3 = D.CrossProd(-C);

        // // Edges cross when the corners have mismatching cross-product
        // // directions.  For a convex shape, opposite corners will be
        // // anti-aligned when the polygon has crossing edges.  A shape with a
        // // concave corner will have three aligned corners, but (I think) no edge
        // // crossings.  So, to find crossing edges, we need to test both pairs of
        // // opposite corners.  If both are anti-aligned then there must be a
        // // crossing edge.
        // if ((corner0.DotProd(corner2) < 0) && (corner1.DotProd(corner3) < 0))
        // {
        //     // If corners 0 and 1 are aligned, then swap 2/3 otherwise swap 1/2
        //     const int swap = (corner0.DotProd(corner1) > 0) ? 3 : 1;
        //     // std::cerr << "\nv0 " << vertices[0]
        //     //           << "\nv1 " << vertices[1]
        //     //           << "\nv2 " << vertices[2]
        //     //           << "\nv3 " << vertices[3]
        //     //           << "\nA " << A
        //     //           << "\nB " << B
        //     //           << "\nC " << C
        //     //           << "\nD " << D
        //     //           << "\nc0 " << corner0
        //     //           << "\nc1 " << corner1
        //     //           << "\nc2 " << corner2
        //     //           << "\nc3 " << corner3
        //     //           << "\nswap 2/" << swap
        //     //           << std::endl;
        //     std::iter_swap(vertices.begin() + 2, vertices.begin() + swap);
        // }
    }
}