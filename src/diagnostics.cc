#include "config.h"

#include <fstream>
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
#include <s2/s2text_format.h>
#include <s2/s2region.h>

#include "constants.h"
#include "ephemeris.h"
#include "util/polygon.h"
#include "util/string.h"
#include "util/spherical.h"

using namespace sbsearch;
using namespace sbsearch::util;
using std::cerr;
using std::endl;

void snapping()
{
    std::cout << "Snapping\n";

    S2Polygon a, b, c, d, e;
    make_polygon(make_vertices("0:0, 1:0, 1:1"), a);
    make_polygon(make_vertices("0:0, 0:1, 1:1"), b);
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
    make_polygon(make_vertices("0:0, 1:0, 1:1"), a);

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

void ellipse_tests()
{
    const int N = 10;
    double ra = 0, dec = 0;
    double mu = 10, mu_theta = 0;
    double unc_a = 180;
    double unc_b = unc_a / 10;
    double step = 0.05;

    Ephemeris eph;

    // for (int i = 0; i < N; i++)
    // {
    //     double unc_theta = (double)std::rand() / RAND_MAX * 360;

    //     eph.append({{(double)i, (double)i,
    //                  ra, dec, mu, mu_theta,
    //                  unc_a, unc_b, unc_theta,
    //                  1, 1, 0, 0, 0, 0, 0, 0}});

    //     ra += step * 1440 * ARCSEC / DEG * mu * std::sin(mu_theta * DEG);
    //     dec += step * 1440 * ARCSEC / DEG * mu * std::cos(mu_theta * DEG);

    //     mu_theta += 360 * step;
    //     mu *= 1 + 1.0 / N;
    // }

    mu_theta = 90;
    for (int i = 0; i < N; i++)
    {
        double unc_theta = (double)std::rand() / RAND_MAX * 360;

        eph.append({{(double)i, (double)i,
                     ra, dec, mu, mu_theta,
                     10 * unc_a, 10 * unc_b, unc_theta,
                     1, 1, 0, 0, 0, 0, 0, 0}});

        ra += step * 1440 * ARCSEC / DEG * mu * std::sin(mu_theta * DEG);
        // dec += step * 1440 * ARCSEC / DEG * mu * std::cos(mu_theta * DEG);

        mu *= 1 - 1.0 / N;
    }

    cerr << eph;
    eph.mutable_options()->use_uncertainty = true;
    auto polygons = eph.as_polygons();

    // std::array<S2Polygon, 2 * N - 1> polygons;
    // for (int i = 0; i < N; i++)
    // {
    //     circumscribe_ellipse(S2LatLng(eph.vertex(i)),
    //                          eph.data(i).unc_a.value() * ARCSEC,
    //                          eph.data(i).unc_b.value() * ARCSEC,
    //                          eph.data(i).unc_theta.value() * DEG,
    //                          eph.data(i).mu_theta.value() * DEG,
    //                          polygons[2 * i]);
    // }

    // // fill in the space between ellipses
    // auto last = polygons.begin();
    // auto next = std::next(polygons.begin(), 2);
    // while (next < polygons.end())
    // {
    //     // Four vertices can build three non-degenerate polygons.  Connect the
    //     // forward edge of the last circumscribing parallelogram to the last
    //     // edge of the next parallelogram.  Check for crossing edges and
    //     // re-order vertices as needed.
    //     std::vector<S2Point> vertices({{last->loop(0)->vertex(0),
    //                                     next->loop(0)->vertex(1),
    //                                     next->loop(0)->vertex(2),
    //                                     last->loop(0)->vertex(3)}});

    //     // Calculate the cross-product of the edges at each corner.
    //     const S2Point A = vertices[1] - vertices[0];
    //     const S2Point B = vertices[2] - vertices[1];
    //     const S2Point C = vertices[3] - vertices[2];
    //     const S2Point D = vertices[0] - vertices[3];
    //     const S2Point corner0 = A.CrossProd(-D);
    //     const S2Point corner1 = B.CrossProd(-A);
    //     const S2Point corner2 = C.CrossProd(-B);
    //     const S2Point corner3 = D.CrossProd(-C);

    //     // Edges cross when the corners have mismatching cross-product
    //     // directions.  For a convex shape, opposite corners will be
    //     // anti-aligned when the polygon has crossing edges.  A shape with a
    //     // concave corner will have three aligned corners, but (I think) no edge
    //     // crossings.  So, to find crossing edges, we need to test both pairs of
    //     // opposite corners.  If both are anti-aligned then there must be a
    //     // crossing edge.
    //     if ((corner0.DotProd(corner2) < 0) && (corner1.DotProd(corner3) < 0))
    //     {
    //         // If corners 0 and 1 are aligned, then swap 2/3 otherwise swap 1/2
    //         int swap = (corner0.DotProd(corner1) > 0) ? 3 : 1;
    //         std::iter_swap(vertices.begin() + 2, vertices.begin() + swap);
    //     }

    //     std::unique_ptr<S2Loop> loop = std::make_unique<S2Loop>(vertices, S2Debug::DISABLE);
    //     loop->Normalize();
    //     std::next(last)->Init(std::move(loop));

    //     std::advance(last, 2);
    //     std::advance(next, 2);
    // }

    std::ofstream outf("polygons.txt", std::ios::binary);

    int start = 0, stop = N;
    cerr << eph.data(0).unc_a.value() * ARCSEC / DEG;
    outf << "ellipses\n";
    for (int i = start; i < stop; i++)
        outf << join<S2LatLng>(
                    ellipse(1000,
                            S2LatLng(eph.vertex(i)),
                            eph.data(i).unc_a.value() * ARCSEC,
                            eph.data(i).unc_b.value() * ARCSEC,
                            eph.data(i).unc_theta.value() * DEG),
                    " ")
             << endl;

    outf << "polygons\n";
    for (int i = 0; i < 2 * stop - 1; i++)
        outf << s2textformat::ToString(*polygons[i]) << endl;

    outf << "ephemeris\n";
    for (int i = start; i < stop; i++)
        outf << eph.data(i).ra.value() << " "
             << eph.data(i).dec.value() << " "
             << eph.data(i).mu_theta.value()
             << endl;

    return;

    for (int i = 0; i < polygons.size(); i++)
    {
        outf << polygons[i]->loop(0)->vertex(0);
    }
}

int main(int argc, char *argv[])
{
    // snapping();
    // buffering();
    ellipse_tests();
    // S2LatLng center = S2LatLng::FromDegrees(0, 0);
    // for (int i = 0; i < 8; i++)
    // {
    //     S1Angle t;
    //     double m;
    //     double a = 1 * DEG;
    //     double b = 0.5 * DEG;
    //     double theta = i * 45 * DEG;
    //     double position_angle = 0;

    //     if (std::abs(position_angle - theta) < 1e-8)
    //     {
    //         t = S1Angle::Degrees(90);
    //     }
    //     else
    //     {
    //         m = std::tan(position_angle - theta);
    //         t = S1Angle::Radians(std::atan2(-b, m * a));
    //     }

    //     // the distance to the tangent points
    //     const S1Angle rho = S1Angle::Radians(ellipse_rho(a, b, t.radians()));

    //     cerr << theta / DEG << " " << t << " " << m << " " << rho << endl;
    // }
    return 0;
}
