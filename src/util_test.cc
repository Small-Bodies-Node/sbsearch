#include "util.h"
#include "sbsearch.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>

#include <s2/s2latlng.h>
#include <s2/s2latlng_rect.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>
#include <s2/s2text_format.h>

#include <gtest/gtest.h>

using std::ceil;
using std::floor;
using std::string;
using std::vector;

namespace sbsearch
{
    TEST(UtilTests, UtilMjdToTimeTerms)
    {
        // test 2 time resolutions, 11 points per time resolution, width of 3 points
        // expect: first 8 steps have one term, next 3 have two terms, repeat
        for (int step = 0; step < 2 * 11; step++)
        {
            const double start = 59800.0 + step / 11.0 / TIME_TERMS_PER_DAY;
            const double stop = 59800.0 + (step + 3) / 11.0 / TIME_TERMS_PER_DAY;
            const vector<string> terms = mjd_to_time_terms(start, stop);

            EXPECT_EQ(terms[0], std::to_string(59800 * TIME_TERMS_PER_DAY + (unsigned int)floor(step / 11.0)));
            if (step % 11 < 9)
            {
                EXPECT_EQ(terms.size(), 1);
            }
            else
            {
                EXPECT_EQ(terms.size(), 2);
                EXPECT_EQ(terms[1], std::to_string(59800 * TIME_TERMS_PER_DAY + (unsigned int)floor((step + 3) / 11.0)));
            }
        }
    }

    TEST(UtilTests, UtilPositionAngle)
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

    TEST(UtilTests, UtilSplit)
    {
        const string s = ",1,22, 3, ";
        const vector<string> parts = split(s, ',');
        const vector<string> expected = {"", "1", "22", " 3", " "};
        EXPECT_EQ(parts, expected);
    }

    TEST(UtilTests, UtilJoin)
    {
        const vector<string> parts = {"", "1", "22", " 3", " "};
        const string s = join(parts, ",");
        const string expected = ",1,22, 3, ";
        EXPECT_EQ(s, expected);
    }

    TEST(UtilTests, UtilFormatVertices)
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

    TEST(UtilTests, UtilMakeVertices)
    {
        string formatted = "0.000000:0.000000, 1.000000:0.000000, 1.000000:1.000000, 0.000000:1.000000";
        vector<S2Point> points = makeVertices(formatted);
        vector<S2Point> expected = {
            S2LatLng::FromDegrees(0, 0).ToPoint(),
            S2LatLng::FromDegrees(0, 1).ToPoint(),
            S2LatLng::FromDegrees(1, 1).ToPoint(),
            S2LatLng::FromDegrees(1, 0).ToPoint()};
        for (int i = 0; i < 4; i++)
            EXPECT_EQ(points[i], expected[i]);
    }

    TEST(UtilTests, UtilMakeVerticesErrors)
    {
        EXPECT_THROW(makeVertices("0, 1:0"), std::runtime_error);
        EXPECT_THROW(makeVertices("0:a, 1:0, 1:1, 0:1"), std::runtime_error);
    }

    TEST(UtilTests, UtilMakePolygon)
    {
        std::unique_ptr<S2Polygon> polygon = makePolygon("0:0, 1:0, 1:1, 0:1");
        // note: s2geometry's text format is lat:lng
        std::unique_ptr<S2Polygon> expected = s2textformat::MakePolygonOrDie("0:0, 0:1, 1:1, 1:0");
        EXPECT_TRUE(polygon->BoundaryEquals(*expected.get()));
    }

    // I'm not sure how to force an error here
    // TEST(UtilTests, UtilMakePolygonErrors)
    // {
    //     EXPECT_THROW(makePolygon("0:0"), std::runtime_error);
    // }
}