// Licensed with the 3-clause BSD license.  See LICENSE for details.

#include <iostream>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <s2/s2builder.h>
#include <s2/s2builderutil_s2polygon_layer.h>
#include <s2/s2cap.h>
#include <s2/s2error.h>
#include <s2/s2loop.h>
#include <s2/s2polygon.h>
#include <s2/s2polyline.h>
#include <s2/s2point.h>

using std::asin;
using std::atan2;
using std::cos;
using std::sin;
using std::tan;
using std::vector;

enum IntersectionType
{
    PolygonContainsPoint = 0,
    PolygonContainsArea = 1,
    PolygonIntersectsArea = 2,
    AreaContainsPolygon = 3,
};

static void _build_polygon(double *, double *, int, bool, S2Polygon &);
static void _build_polygon_from_vertices(vector<S2Point>, bool, S2Polygon &);
static S2Polyline _build_polyline(double *, double *, int, double, double);
static bool _verify_build_polygon(double *, double *, int);
static bool _polygon_contains_point(double *, double *, int, double, double);
static bool _polygon_intersects_line(double *, double *, int, double *, double *, int, double, double);
static bool _polygon_intersects_about_line(double *, double *, int, double *, double *, int, double *, double *, double, double);
static bool _polygon_intersects_polygon(double *, double *, int, double *, double *, int);
static bool _polygon_intersects_cap(double *, double *, int, double, double, double, int);

// Meeus 1998, Astronomical Algorithms, p116
static double _position_angle(double ra1, double dec1, double ra2, double dec2)
{
    double dra = ra2 - ra1;
    return atan2(sin(dra), cos(dec1) * tan(dec2) - sin(dec1) * cos(dra));
}

// Algorithm based on astropy.coordinates.angle_utilities.offset_by
static void _offset_by(double ra, double dec, double pa, double rho, double &new_ra, double &new_dec)
{
    double cos_a = cos(rho);
    double sin_a = sin(rho);
    // note: c is measured from the pole, dec from equator
    double cos_c = sin(dec);
    double sin_c = cos(dec);
    double cos_B = cos(pa);
    double sin_B = sin(pa);

    double cos_b = cos_c * cos_a + sin_c * sin_a * cos_B;
    double xsin_A = sin_a * sin_B * sin_c;
    double xcos_A = cos_a - cos_b * cos_c;

    double A = atan2(xsin_A, xcos_A);

    // handle the pole as a very small angle?
    if (sin_c < 1e-12)
        A = M_PI_2 + cos_c * (M_PI_2 - pa);

    new_ra = ra + A;
    new_dec = asin(cos_b);
}

static S2Point offset_point_by(S2Point point, double pa, double rho)
{
    double new_ra, new_dec;
    S2LatLng ll(point);
    _offset_by(ll.coords()[1], ll.coords()[0], pa, rho, new_ra, new_dec);
    return S2LatLng::FromRadians(new_dec, new_ra).Normalized().ToPoint();
}

/* Build polygon from vertices.

Uses S2Builder in order to accomodate loops.

The vertices must form a closed shape, e.g., last vertex = first vertex.
Use close=True if they do not.
*/
static void _build_polygon(double *ra, double *dec, int n, bool close, S2Polygon &polygon)
{
    vector<S2Point> vertices;
    for (int i = 0; i < n; i++)
        vertices.push_back(S2LatLng::FromRadians(dec[i], ra[i]).Normalized().ToPoint());

    _build_polygon_from_vertices(vertices, close, polygon);
}

static void _build_polygon_from_vertices(vector<S2Point> vertices, bool close, S2Polygon &polygon)
{
    S2Builder::Options builder_options;
    builder_options.set_split_crossing_edges(true);
    S2Builder builder(builder_options);

    s2builderutil::S2PolygonLayer::Options layer_options;
    layer_options.set_edge_type(S2Builder::EdgeType::UNDIRECTED);
    builder.StartLayer(std::make_unique<s2builderutil::S2PolygonLayer>(&polygon, layer_options));

    int n = vertices.size();
    for (int i = 0; i < n - 1; i++)
        builder.AddEdge(vertices[i], vertices[i + 1]);

    if (close)
        builder.AddEdge(vertices.back(), vertices.front());

    S2Error error;
    builder.Build(&error);
    if (!error.ok())
    {
        std::string msg = std::string("PolygonBuildError (") + std::to_string(error.code()) + std::string("): ") + std::string(error.text());
        throw std::invalid_argument(msg);
    }
}

static S2Polyline _build_polyline(double *ra, double *dec, int n, double line_start = 0, double line_stop = 1)
{
    vector<S2Point> vertices;
    for (int i = 0; i < n; i++)
        vertices.push_back(S2LatLng::FromRadians(dec[i], ra[i]).Normalized().ToPoint());
    S2Polyline line(vertices);

    if ((line_start != 0) & (line_stop != 1))
    {
        // interpolate to a sub-set
        S2Point start = line.Interpolate(line_start);
        S2Point stop = line.Interpolate(line_stop);
        int next_vertex, last_vertex;
        vector<S2Point> new_vertices;

        line.GetSuffix(line_start, &next_vertex);
        line.GetSuffix(line_stop, &last_vertex);

        // build the new line
        new_vertices.push_back(start);
        if (next_vertex != last_vertex)
            for (int i = last_vertex; i < next_vertex; i++)
                new_vertices.push_back(vertices[i]);
        new_vertices.push_back(stop);

        line = S2Polyline(new_vertices);
    }

    return line;
}

/* For testing.

   Does not close the polygon.
*/
static bool _verify_build_polygon(double *ra, double *dec, int n)
{
    S2Polygon polygon;
    _build_polygon(ra, dec, n, false, polygon);
    return polygon.IsValid();
}

/* Test if the polygon covers the point.

    Parameters
    ----------
    poly_ra, poly_dec : ndarray
        Polygon RA and Dec, radians.

    point_ra, point_dec : ndarray
        Point RA and Dec, radians.


    Returns
    -------
    intersects : bool

*/
static bool _polygon_contains_point(double *poly_ra, double *poly_dec, int n,
                                   double point_ra, double point_dec)
{

    S2Polygon polygon;
    _build_polygon(poly_ra, poly_dec, n, true, polygon);

    S2Point point = S2LatLng::FromRadians(point_dec, point_ra).Normalized().ToPoint();

    return polygon.Contains(point);
}

/*  Test if the polygon intersects the line.


    Parameters
    ----------
    poly_ra, poly_dec : ndarray
        Polygon RA and Dec, radians.

    line_ra, line_dec : ndarray
        Line RA and Dec, radians.

    line_start, line_stop : float
        Test a sub-portion of the line, indicated by the given line fractions
        (linear interpolation).


    Returns
    -------
    intersects : bool

*/
static bool _polygon_intersects_line(double *poly_ra, double *poly_dec, int poly_n,
                                    double *line_ra, double *line_dec, int line_n,
                                    double line_start = 0, double line_stop = 1)
{
    if (line_stop <= line_start)
        throw std::domain_error("line_stop <= line_start");

    S2Polygon polygon;
    _build_polygon(poly_ra, poly_dec, poly_n, true, polygon);

    S2Polyline line = _build_polyline(line_ra, line_dec, line_n, line_start, line_stop);

    return polygon.Intersects(line);
}

/* Test for intersection between polygon and region about a line.


    Parameters
    ----------
    poly_ra, poly_dec : ndarray
        Polygon RA and Dec, radians.

    line_ra, line_dec : ndarray
        Line RA and Dec, radians.

    line_start, line_stop : float
        Test a sub-portion of the line, indicated by the given line fractions.

    a, b : ndarray
        Angular size of the region: ``a`` is along the line, ``b``
        is perpendicular to it, radians.  Only the first and last
        ``a`` are considered.

    line_start, line_stop : float
        Test a sub-portion of the line, indicated by the given line fractions
        (linear interpolation).


    Returns
    -------
    intersects : bool

*/
static bool _polygon_intersects_about_line(double *poly_ra, double *poly_dec, int poly_n,
                                          double *line_ra, double *line_dec, int line_n,
                                          double *a, double *b,
                                          double line_start = 0, double line_stop = 1)
{
    S2Polyline line = _build_polyline(line_ra, line_dec, line_n, line_start, line_stop);

    // transform the line into a polygon
    // the spine is the input line extended by + / -a
    // the polygon is the region betweeh spine + b and spine - b
    vector<double> pa;
    vector<S2Point> spine;

    for (int i = 0; i < line_n - 1; i++)
        pa.insert(pa.end(), _position_angle(line_ra[i], line_dec[i], line_ra[i + 1], line_dec[i + 1]));
    pa.insert(pa.begin(), pa.front());
    pa.insert(pa.end(), pa.back());
    pa.insert(pa.end(), pa.back());

    for (int i = 0; i < line_n; i++)
        spine.insert(spine.end(), line.vertex(i));
    spine.insert(spine.begin(), offset_point_by(line.vertex(0), pa.front(), -a[0]));
    spine.insert(spine.end(), offset_point_by(line.vertex(line.num_vertices() - 1), pa.back(), a[line_n - 1]));

    // Construct query region polygon. In S2, the interior of a line is on the
    // left.  This is CCW from each edge for small loops.  S2 is designed for
    // the surface of the Earth, but equatorial coordinates on the sky are a
    // mirror of the Earth's coordinates.  Therefore, the interior is on the
    // right (the CW edge).  Thus, we first generate the edges along PA + pi /
    // 2, then return to the first vertex along PA - pi / 2.
    vector<S2Point> vertices;
    for (int i = 0; i < line_n + 2; i++)
    {
        // index for pa, a, and b
        int j = i - 1;
        j = static_cast<int>(std::max(0, j));
        j = static_cast<int>(std::min(j, line_n - 1));
        vertices.insert(vertices.end(), offset_point_by(spine[i], pa[i] + M_PI_2, b[j]));
        vertices.insert(vertices.begin(), offset_point_by(spine[i], pa[i] - M_PI_2, b[j]));
    }

    S2Polygon polygon1, polygon2;
    _build_polygon(poly_ra, poly_dec, line_n + 2, true, polygon1);
    _build_polygon_from_vertices(vertices, true, polygon2);

    return polygon1.Intersects(&polygon2);
}

// Test for polygon intersection.
static bool _polygon_intersects_polygon(double *ra1, double *dec1, int n1,
                                        double *ra2, double *dec2, int n2)
{
    S2Polygon polygon1, polygon2;
    _build_polygon(ra1, dec1, n1, true, polygon1);
    _build_polygon(ra2, dec2, n2, true, polygon2);
    return polygon1.Intersects(polygon2);
}

static bool _polygon_intersects_cap(double *poly_ra, double *poly_dec, int poly_n,
                                   double point_ra, double point_dec, double radius,
                                   int intersection)
{
    S2Polygon polygon;
    _build_polygon(poly_ra, poly_dec, poly_n, true, polygon);

    S2Cap cap(S2LatLng::FromRadians(point_dec, point_ra).Normalized().ToPoint(), S1Angle::Radians(radius));

    IntersectionType itype = static_cast<IntersectionType>(intersection);

    bool result = false;

    switch (intersection)
    {
    case PolygonContainsPoint:
        result = polygon.Contains(cap.center());
        break;
    case PolygonContainsArea:
        result = (
            polygon.GetDistanceToBoundary(cap.center()) > cap.radius().ToAngle()
            & polygon.Contains(cap.center())
        );
        break;
    case PolygonIntersectsArea:
        result = polygon.GetDistance(cap.center()) < cap.radius().ToAngle();
        break;
    case AreaContainsPolygon:
        // only testing loop[0]; sbsearch does not support multiple loops
        S2Loop *loop = polygon.loop(0);
        result = true;
        for (int i = 0; i < loop->num_vertices(); i++)
        {
            if (!cap.InteriorContains(loop->vertex(i)))
            {
                result = false;
                break;
            }
        }
        break;
    }

    return result;
}