#include "config.h"

#include <stdexcept>
#include <string>

#include <gtest/gtest.h>
#include <s2/s2latlng.h>
#include <s2/s2latlng_rect.h>
#include <s2/s2point.h>

#include "string.h"

using std::string;
using std::vector;

namespace sbsearch::util
{
    TEST(UtilStringTests, Split)
    {
        vector<string> parts = split(",1,22, 3, ", ',');
        vector<string> expected = {"", "1", "22", " 3", " "};
        EXPECT_EQ(parts, expected);

        parts = split("a,b,casdf", ',');
        expected = {"a", "b", "casdf"};
        EXPECT_EQ(parts, expected);

        parts = split("a,b,casdf,", ','); // delimiter terminated
        expected = {"a", "b", "casdf"};
        EXPECT_EQ(parts, expected);
    }

    TEST(UtilStringTests, JoinString)
    {
        const vector<string> parts = {"", "1", "22", " 3", " "};
        const string s = join(parts, ",");
        EXPECT_EQ(s, ",1,22, 3, ");
    }

    TEST(UtilStringTests, JoinDouble)
    {
        const vector<double> parts = {1, 2, 3, 55.5};
        const string s = join(parts, ",");
        EXPECT_EQ(s, "1,2,3,55.5");
    }

    TEST(UtilStringTests, FormatVertices)
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

        EXPECT_EQ(format_vertices(latlngs), formatted);
        EXPECT_EQ(format_vertices(points), formatted);
        EXPECT_EQ(format_vertices(rect), formatted);
    }

    TEST(UtilStringTests, MakeVertices)
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

    TEST(UtilStringTests, MakeVerticesErrors)
    {
        EXPECT_THROW(make_vertices("0, 1:0"), std::runtime_error);
        EXPECT_THROW(make_vertices("0:a, 1:0, 1:1, 0:1"), std::runtime_error);
    }
}
