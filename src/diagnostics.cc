#include "config.h"

#include <iostream>
#include <memory>
#include <s2/s1angle.h>
#include <s2/s1chord_angle.h>
#include <s2/s2buffer_operation.h>
#include <s2/s2builder.h>
#include <s2/s2builderutil_lax_polygon_layer.h>
#include <s2/s2builderutil_snap_functions.h>
#include <s2/s2cap.h>
#include <s2/s2latlng.h>
#include <s2/s2lax_polygon_shape.h>
#include <s2/s2polygon.h>
#include <s2/s2region.h>

#include "util.h"

using namespace sbsearch;

void snapping()
{
    std::cout << "Snapping\n";

    S2Polygon a, b, c, d, e;
    make_polygon("0:0, 1:0, 1:1", a);
    make_polygon("0:0, 0:1, 1:1", b);
    c.InitToIntersection(a, b);

    S1Angle tolerance = S1Angle::Degrees(0.2);
    d.InitToIntersection(a, b, s2builderutil::IdentitySnapFunction(tolerance));

    std::cout << "a, b, a∩b, a∩b with 0.2 deg snap:\n"
              << format_vertices(a) << "\n"
              << format_vertices(b) << "\n"
              << format_vertices(c) << "\n"
              << format_vertices(d) << "\n";
}

void buffering()
{
    std::cout << "Buffering\n";

    S2Polygon a, b, c, d, e;
    make_polygon("0:0, 1:0, 1:1", a);

    S2BufferOperation::Options buffer_options(S1Angle::Degrees(0.2));

    auto output = std::make_unique<S2LaxPolygonShape>();
    S2BufferOperation op(std::make_unique<s2builderutil::LaxPolygonLayer>(output.get()),
                         buffer_options);
    S2Polygon::Shape shape(&a);
    op.AddShape(shape);

    S2Error error;
    if (!op.Build(&error))
    {
        std::cerr << error.text() << std::endl;
        return;
    }

    vector<S2Point> vertices;
    for (int i = 0; i < output->num_loop_vertices(0); i++)
        vertices.push_back(output->loop_vertex(0, i));
    make_polygon(vertices, b);

    // shrink
    buffer_options.set_buffer_radius(S1Angle::Degrees(-0.2));
    auto output2 = std::make_unique<S2LaxPolygonShape>();
    S2BufferOperation op2(std::make_unique<s2builderutil::LaxPolygonLayer>(output2.get()),
                          buffer_options);
    op2.AddShape(shape);

    if (!op2.Build(&error))
    {
        std::cerr << error.text() << std::endl;
        return;
    }

    vertices.clear();
    for (int i = 0; i < output2->num_loop_vertices(0); i++)
        vertices.push_back(output2->loop_vertex(0, i));
    make_polygon(vertices, c);

    // large buffer
    buffer_options.set_buffer_radius(S1Angle::Degrees(2));
    auto output4 = std::make_unique<S2LaxPolygonShape>();
    S2BufferOperation op4(std::make_unique<s2builderutil::LaxPolygonLayer>(output4.get()),
                          buffer_options);
    op4.AddShape(shape);

    if (!op4.Build(&error))
    {
        std::cerr << error.text() << std::endl;
        return;
    }

    vertices.clear();
    for (int i = 0; i < output4->num_loop_vertices(0); i++)
        vertices.push_back(output4->loop_vertex(0, i));
    make_polygon(vertices, d);

    // repeat buffering, but with the lowest resolution allowed
    buffer_options.set_buffer_radius(S1Angle::Degrees(0.2));
    buffer_options.set_circle_segments(3);
    auto output3 = std::make_unique<S2LaxPolygonShape>();
    S2BufferOperation op3(std::make_unique<s2builderutil::LaxPolygonLayer>(output3.get()),
                          buffer_options);
    op3.AddShape(shape);

    if (!op3.Build(&error))
    {
        std::cerr << error.text() << std::endl;
        return;
    }

    vertices.clear();
    for (int i = 0; i < output3->num_loop_vertices(0); i++)
        vertices.push_back(output3->loop_vertex(0, i));
    make_polygon(vertices, e);

    std::cout << "a, a with 0.2 deg buffer with defaults, then shrink by the same amount, then buffer by a large amount, then buffering with the lowest resolution allowed:\n"
              << format_vertices(a) << ";"
              << format_vertices(b) << ";"
              << format_vertices(c) << ";"
              << format_vertices(d) << ";"
              << format_vertices(e) << "\n";
}

int main(int argc, char *argv[])
{
    snapping();
    buffering();
    return 0;
}
