#include "config.h"

#include <stdexcept>
#include <string>
#include <gtest/gtest.h>
#include <s2/s2latlng.h>
#include <s2/s2latlng_rect.h>
#include <s2/s2point.h>

#include "vertices.h"

using std::string;
using std::vector;

namespace sbsearch::util
{
    TEST(VerticesTests, FormatVertices)
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

    TEST(VerticesTests, MakeVertices)
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

    TEST(VerticesTests, MakeVerticesErrors)
    {
        EXPECT_THROW(make_vertices("0, 1:0"), std::runtime_error);
        EXPECT_THROW(make_vertices("0:a, 1:0, 1:1, 0:1"), std::runtime_error);
    }
}
