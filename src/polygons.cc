#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>
#include <s2/s1angle.h>
#include <s2/s2builder.h>
#include <s2/s2buffer_operation.h>
#include <s2/s2builder.h>
#include <s2/s2builderutil_lax_polygon_layer.h>
#include <s2/s2builderutil_s2polygon_layer.h>
#include <s2/s2debug.h>
#include <s2/s2error.h>
#include <s2/s2latlng.h>
#include <s2/s2lax_polygon_shape.h>
#include <s2/s2loop.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>

#include "config.h"
#include "constants.h"
#include "ephemeris/ephemeris.h"
#include "logging.h"
#include "polygons.h"
#include "util/math.h"
#include "util/spherical.h"
#include "util/string.h"

using sbsearch::ephemeris::Ephemeris;

namespace sbsearch
{
    void make_polygon_simple(const vector<S2Point> &vertices, S2Polygon &polygon)
    {
        polygon.Release();
        std::unique_ptr<S2Loop> loop = std::make_unique<S2Loop>(vertices,
                                                                S2Debug::DISABLE);
        loop->Normalize();
        polygon.Init(std::move(loop));
    }

    void make_polygon(const vector<S2Point> &vertices, S2Polygon &polygon)
    {
        polygon.Release();

        S2Builder::Options builder_options;
        builder_options.set_split_crossing_edges(true);
        S2Builder builder{builder_options};

        s2builderutil::S2PolygonLayer::Options layer_options;
        layer_options.set_edge_type(S2Builder::EdgeType::UNDIRECTED);
        builder.StartLayer(
            std::make_unique<s2builderutil::S2PolygonLayer>(&polygon, layer_options));

        int n = vertices.size();
        for (int i = 0; i < n - 1; i++)
            builder.AddEdge(vertices[i], vertices[i + 1]);

        // close the polygon
        builder.AddEdge(vertices[n - 1], vertices[0]);

        S2Error error;
        builder.Build(&error);
        if (!error.ok())
        {
            Logger::error() << error.code() << " " << error.text() << std::endl;
            throw std::runtime_error("Polygon build error: " + error.text());
        }
    }

    void padded_polygon(const S2Polygon &polygon,
                        const double padding,
                        S2Polygon &result)
    {
        if (padding <= 0)
        {
            // nothing to do
            result.Copy(polygon);
            return;
        }

        S2BufferOperation::Options buffer_options(S1Angle::Degrees(padding / 60));

        auto output = std::make_unique<S2LaxPolygonShape>();
        S2BufferOperation op(
            std::make_unique<s2builderutil::LaxPolygonLayer>(output.get()),
            buffer_options);
        S2Polygon::Shape shape(&polygon);
        op.AddShape(shape);

        S2Error error;
        if (!op.Build(&error))
        {
            Logger::error() << error.code() << " " << error.text() << std::endl;
            throw std::runtime_error("Polygon buffer error: " + error.text());
        }

        vector<S2Point> vertices;
        for (int i = 0; i < output->num_loop_vertices(0); i++)
            vertices.push_back(output->loop_vertex(0, i));

        make_polygon(vertices, result);
        return;
    }

    void make_polygon(const Observation &observation,
                      S2Polygon &polygon,
                      const bool verify)
    {
        polygon.set_s2debug_override(S2Debug::DISABLE);
        if (verify)
            make_polygon(util::make_vertices(observation.fov()), polygon);
        else
            make_polygon_simple(util::make_vertices(observation.fov()), polygon);
    }

    vector<unique_ptr<S2Polygon>> make_polygons(const Ephemeris &eph,
                                                const bool use_uncertainty,
                                                double padding)
    {
        // The minimum padding is a 0.1" radius circle
        padding = std::max(padding * 60.0, 0.1);
        vector<double> a(eph.num_vertices(), padding),
            b(eph.num_vertices(), padding),
            theta(eph.num_vertices(), 0);

        if (use_uncertainty)
        {
            for (int i = 0; i < eph.num_vertices(); i++)
            {
                a[i] += eph.data()[i].unc_a.value_or(0);
                b[i] += eph.data()[i].unc_b.value_or(0);
                theta[i] = eph.data()[i].unc_theta.value_or(0);
            }
        }

        // collect polygons in a vector
        std::vector<std::unique_ptr<S2Polygon>> polygons;
        polygons.reserve(eph.num_vertices() + eph.num_segments());

        // polygons around each ephemeris point
        for (int i = 0; i < eph.num_vertices(); i++)
            polygons.emplace_back(
                util::circumscribe_ellipse(
                    S2LatLng(eph.vertex(i)),
                    a[i] * ARCSEC,
                    b[i] * ARCSEC,
                    theta[i] * DEG,
                    eph.data()[i].mu_theta.value() * DEG));

        // polygons between ephemeris points
        std::vector<S2Point> vertices(4);
        for (int i = 0; i < eph.num_vertices() - 1; i++)
        {
            vertices = {{polygons[i]->loop(0)->vertex(0),
                         polygons[i + 1]->loop(0)->vertex(1),
                         polygons[i + 1]->loop(0)->vertex(2),
                         polygons[i]->loop(0)->vertex(3)}};
            util::fix_crossing_edges(vertices);
            auto loop = std::make_unique<S2Loop>(vertices, S2Debug::DISABLE);
            loop->Normalize();
            auto polygon = std::make_unique<S2Polygon>(S2Polygon(std::move(loop),
                                                                 S2Debug::DISABLE));
            polygons.emplace_back(std::move(polygon));
        }

        return std::move(polygons);
    }
}