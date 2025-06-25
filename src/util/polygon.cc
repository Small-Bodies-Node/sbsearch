#include "config.h"

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

#include "polygon.h"
#include "../logging.h"

using std::vector;

namespace sbsearch::util
{
    void make_polygon_simple(const vector<S2Point> &vertices, S2Polygon &polygon)
    {
        polygon.Release();
        std::unique_ptr<S2Loop> loop = std::make_unique<S2Loop>(vertices, S2Debug::DISABLE);
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
        builder.StartLayer(std::make_unique<s2builderutil::S2PolygonLayer>(&polygon, layer_options));

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

    void padded_polygon(const S2Polygon &polygon, const double padding, S2Polygon &result)
    {
        if (padding <= 0)
        {
            // nothing to do
            result.Copy(polygon);
            return;
        }

        S2BufferOperation::Options buffer_options(S1Angle::Degrees(padding / 60));

        auto output = std::make_unique<S2LaxPolygonShape>();
        S2BufferOperation op(std::make_unique<s2builderutil::LaxPolygonLayer>(output.get()),
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
}