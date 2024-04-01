#include "config.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>
#include <s2/s2latlng.h>
#include <s2/s2latlng_rect.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>
#include <s2/s2text_format.h>

#include "util.h"

using std::ceil;
using std::floor;
using std::string;
using std::vector;

namespace sbsearch
{
    TEST(UtilTests, PositionAngle)
    {
        /*
        # verify with astropy
        >>> a = SkyCoord('0d 0d')
        >>> c = SkyCoord('0.01d 0.01d')
        >>> a.position_angle(c).deg
        44.999999563667686
        */
        const S2Point a = S2LatLng::FromDegrees(0, 0).ToPoint();
        const S2Point b = S2LatLng::FromDegrees(0, 0.01).ToPoint();
        const S2Point c = S2LatLng::FromDegrees(0.01, 0.01).ToPoint();
        const S2Point d = S2LatLng::FromDegrees(0.01, 0).ToPoint();
        const S2Point e = S2LatLng::FromDegrees(0, -0.01).ToPoint();
        const S2Point f = S2LatLng::FromDegrees(-0.01, 0).ToPoint();
        EXPECT_EQ(position_angle(a, b), 90 * DEG);
        EXPECT_EQ(position_angle(a, c), 44.999999563667686 * DEG);
        EXPECT_EQ(position_angle(a, d), 0);
        EXPECT_EQ(position_angle(a, e), -90 * DEG);
        EXPECT_EQ(position_angle(a, f), 180 * DEG);
    }

    TEST(UtilTests, OffsetBy)
    {
        /*
        # verify with astropy
        >>> SkyCoord(0, 0, unit='rad').directional_offset_by(30 * u.deg, 2 * u.deg)
        <SkyCoord (ICRS): (ra, dec) in deg
            (1.00030471, 1.73196284)>
        >>> SkyCoord(0, 0, unit='rad').directional_offset_by(30 * u.deg, 2 * u.deg).ra.value
        1.0003047102322884
        >>> SkyCoord(0, 0, unit='rad').directional_offset_by(30 * u.deg, 2 * u.deg).dec.value
        1.7319628412802042
        */
        const S2LatLng coords = S2LatLng::FromDegrees(0, 0);

        const S2LatLng a = offset_by(coords, S1Angle::Degrees(0), S1Angle::Degrees(10));
        EXPECT_NEAR(a.lat().degrees(), 10, 1e-8);
        EXPECT_NEAR(a.lng().degrees(), 0, 1e-8);

        const S2LatLng b = offset_by(coords, S1Angle::Degrees(90), S1Angle::Degrees(10));
        EXPECT_NEAR(b.lat().degrees(), 0, 1e-8);
        EXPECT_NEAR(b.lng().degrees(), 10, 1e-8);

        const S2LatLng c = offset_by(coords, S1Angle::Degrees(180), S1Angle::Degrees(10));
        EXPECT_NEAR(c.lat().degrees(), -10, 1e-8);
        EXPECT_NEAR(c.lng().degrees(), 0, 1e-8);

        const S2LatLng d = offset_by(coords, S1Angle::Degrees(270), S1Angle::Degrees(10));
        EXPECT_NEAR(d.lat().degrees(), 0, 1e-8);
        EXPECT_NEAR(d.lng().degrees(), -10, 1e-8);

        const S2LatLng e = offset_by(coords, S1Angle::Degrees(-90), S1Angle::Degrees(10));
        EXPECT_NEAR(e.lat().degrees(), 0, 1e-8);
        EXPECT_NEAR(e.lng().degrees(), -10, 1e-8);

        const S2LatLng f = offset_by(coords, S1Angle::Degrees(30), S1Angle::Degrees(2));
        EXPECT_NEAR(f.lat().degrees(), 1.7319628412802042, 1e-8);
        EXPECT_NEAR(f.lng().degrees(), 1.0003047102322884, 1e-8);
    }

    TEST(UtilTests, Ellipse)
    {
        vector<S2LatLng> e = ellipse(4, S2LatLng::FromRadians(0, 0), 0.1, 0.05, 0);
        EXPECT_NEAR(e[0].lat().radians(), 0.1, 1e-8);
        EXPECT_NEAR(e[0].lng().radians(), 0, 1e-8);
        EXPECT_NEAR(e[1].lat().radians(), 0, 1e-8);
        EXPECT_NEAR(e[1].lng().radians(), 0.05, 1e-8);
        EXPECT_NEAR(e[2].lat().radians(), -0.1, 1e-8);
        EXPECT_NEAR(e[2].lng().radians(), 0, 1e-8);
        EXPECT_NEAR(e[3].lat().radians(), 0, 1e-8);
        EXPECT_NEAR(e[3].lng().radians(), -0.05, 1e-8);

        e = ellipse(4, S2LatLng::FromRadians(0, 0), 0.1, 0.05, PI_2);
        EXPECT_NEAR(e[0].lat().radians(), 0, 1e-8);
        EXPECT_NEAR(e[0].lng().radians(), 0.1, 1e-8);
        EXPECT_NEAR(e[1].lat().radians(), -0.05, 1e-8);
        EXPECT_NEAR(e[1].lng().radians(), 0, 1e-8);
        EXPECT_NEAR(e[2].lat().radians(), 0, 1e-8);
        EXPECT_NEAR(e[2].lng().radians(), -0.1, 1e-8);
        EXPECT_NEAR(e[3].lat().radians(), 0.05, 1e-8);
        EXPECT_NEAR(e[3].lng().radians(), 0, 1e-8);
    }

    TEST(UtilTests, Split)
    {
        const string s = ",1,22, 3, ";
        const vector<string> parts = split(s, ',');
        const vector<string> expected = {"", "1", "22", " 3", " "};
        EXPECT_EQ(parts, expected);
    }

    TEST(UtilTests, Join)
    {
        const vector<string> parts = {"", "1", "22", " 3", " "};
        const string s = join(parts, ",");
        const string expected = ",1,22, 3, ";
        EXPECT_EQ(s, expected);
    }

    TEST(UtilTests, FormatVertices)
    {
        string formatted = "0.000000:0.000000, 1.000000:0.000000, 1.000000:1.000000, 0.000000:1.000000";
        vector<S2LatLng> latlngs = {
            S2LatLng::FromDegrees(0, 0),
            S2LatLng::FromDegrees(0, 1),
            S2LatLng::FromDegrees(1, 1),
            S2LatLng::FromDegrees(1, 0)};
        vector<S2Point> points = {
            S2LatLng::FromDegrees(0, 0).ToPoint(),
            S2LatLng::FromDegrees(0, 1).ToPoint(),
            S2LatLng::FromDegrees(1, 1).ToPoint(),
            S2LatLng::FromDegrees(1, 0).ToPoint()};
        S2LatLngRect rect = S2LatLngRect::FromCenterSize(S2LatLng::FromDegrees(0.5, 0.5), S2LatLng::FromDegrees(1, 1));
        double ra[4] = {0, 1, 1, 0};
        double dec[4] = {0, 0, 1, 1};

        EXPECT_EQ(format_vertices(latlngs), formatted);
        EXPECT_EQ(format_vertices(points), formatted);
        EXPECT_EQ(format_vertices(rect), formatted);
        EXPECT_EQ(format_vertices(4, ra, dec), formatted);
    }

    TEST(UtilTests, MakeVertices)
    {
        string formatted = "0.000000:0.000000, 1.000000:0.000000, 1.000000:1.000000, 0.000000:1.000000";
        vector<S2Point> points = make_vertices(formatted);
        vector<S2Point> expected = {
            S2LatLng::FromDegrees(0, 0).ToPoint(),
            S2LatLng::FromDegrees(0, 1).ToPoint(),
            S2LatLng::FromDegrees(1, 1).ToPoint(),
            S2LatLng::FromDegrees(1, 0).ToPoint()};
        for (int i = 0; i < 4; i++)
            EXPECT_EQ(points[i], expected[i]);
    }

    TEST(UtilTests, MakeVerticesErrors)
    {
        EXPECT_THROW(make_vertices("0, 1:0"), std::runtime_error);
        EXPECT_THROW(make_vertices("0:a, 1:0, 1:1, 0:1"), std::runtime_error);
    }

    TEST(UtilTests, MakePolygon)
    {
        S2Polygon polygon;
        make_polygon("0:0, 1:0, 1:1, 0:1", polygon);
        // note: s2geometry's text format is lat:lng
        std::unique_ptr<S2Polygon> expected = s2textformat::MakePolygonOrDie("0:0, 0:1, 1:1, 1:0");
        EXPECT_TRUE(polygon.BoundaryEquals(*expected.get()));
    }

    // I'm not sure how to force an error here
    // TEST(UtilTests, MakePolygonErrors)
    // {
    //     EXPECT_THROW(makePolygon("0:0"), std::runtime_error);
    // }

    TEST(UtilTests, PaddedPolygon)
    {
        // results based on diagonstics.cc

        S2Polygon polygon, result;
        make_polygon("0:0, 1:0, 1:1", polygon);
        padded_polygon(polygon, 0.2 * 60, result);
        EXPECT_EQ(format_vertices(result), "0.000000:-0.200000, 1.000000:-0.200000, 1.078008:-0.184160, 1.143658:-0.139148, 1.186553:-0.072096, 1.200000:0.000000, 1.200030:0.999994, 1.184192:1.078002, 1.139176:1.143655, 1.072111:1.186553, 0.993622:1.199898, 0.916143:1.181578, 0.858561:1.141429, -0.141432:0.141410, -0.185386:0.075047, -0.199974:-0.003204, -0.182887:-0.080947, -0.136830:-0.145869, -0.069098:-0.187684");

        padded_polygon(polygon, 2.0 * 60, result);
        EXPECT_EQ(format_vertices(result), "0.000000:-2.000000, 1.000000:-2.000000, 1.780631:-1.841419, 2.437301:-1.390889, 2.865971:-0.719955, 3.000000:0.000000, 3.000304:0.999391, 2.842308:1.779711, 2.392148:2.436582, 1.720856:2.865745, 0.934903:2.998942, 0.159309:2.814960, -0.415218:2.413873, -1.414608:1.413962, -1.854083:0.750049, -1.999733:-0.032661, -1.828606:-0.810198, -1.367674:-1.459408, -0.689756:-1.877340");
    }

    TEST(UtilTests, Mjd2Cal)
    {
        EXPECT_EQ(mjd2cal(60000.0), "2023-02-25 00:00");
    }
}
