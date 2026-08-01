#include <memory>
#include <vector>
#include <s2/s2loop.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>

#include "config.h"
#include "constants.h"
#include "ephemeris.h"
#include "util/math.h"
#include "util/spherical.h"

namespace sbsearch
{
    vector<unique_ptr<S2Polygon>> Ephemeris::as_polygons(double padding) const
    {
        // The minimum padding is a 0.1" radius circle
        padding = std::max(padding * 60.0, 0.1);
        vector<double> a(num_vertices_, padding), b(num_vertices_, padding), theta(num_vertices_, 0);

        if (options_.use_uncertainty)
        {
            for (int i = 0; i < num_vertices_; i++)
            {
                a[i] += data_[i].unc_a.value_or(0);
                b[i] += data_[i].unc_b.value_or(0);
                theta[i] = data_[i].unc_theta.value_or(0);
            }
        }

        // collect polygons in a vector
        std::vector<std::unique_ptr<S2Polygon>> polygons;

        // polygons around each ephemeris point
        for (int i = 0; i < num_vertices_; i++)
            polygons.emplace_back(
                util::circumscribe_ellipse(
                    S2LatLng(vertex(i)),
                    a[i] * ARCSEC,
                    b[i] * ARCSEC,
                    theta[i] * DEG,
                    data_[i].mu_theta.value() * DEG));

        // polygons between ephemeris points
        for (int i = 0; i < num_vertices_ - 1; i++)
        {
            std::vector<S2Point> vertices({{polygons[i]->loop(0)->vertex(0),
                                            polygons[i + 1]->loop(0)->vertex(1),
                                            polygons[i + 1]->loop(0)->vertex(2),
                                            polygons[i]->loop(0)->vertex(3)}});
            util::fix_crossing_edges(vertices);
            std::unique_ptr<S2Loop> loop = std::make_unique<S2Loop>(vertices, S2Debug::DISABLE);
            loop->Normalize();
            polygons.emplace_back(std::move(std::make_unique<S2Polygon>(S2Polygon(std::move(loop)))));
        }

        return std::move(polygons);
    }
}