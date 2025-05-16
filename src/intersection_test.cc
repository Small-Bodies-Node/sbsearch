#include <gtest/gtest.h>
#include <s2/s1chord_angle.h>
#include <s2/s2cap.h>
#include <s2/s2latlng.h>
#include <s2/s2polygon.h>

#include "intersection.h"
#include "util/polygon.h"
#include "util/string.h"

using std::vector;

namespace sbsearch
{
    namespace testing
    {
        TEST(SBSearchTests, PolygonIntersectsCap)
        {
            S2Polygon polygon;
            S2Cap area;

            util::make_polygon(util::make_vertices("0:0, 0:1, 1:1, 1:0"), polygon);

            // does not intersect
            area = S2Cap(S2LatLng::FromDegrees(1.1, 1.1).Normalized().ToPoint(), S1ChordAngle::Degrees(0.05));
            EXPECT_FALSE(intersects(polygon, area, ContainsCenter));
            EXPECT_FALSE(intersects(polygon, area, IntersectsArea));
            EXPECT_FALSE(intersects(polygon, area, ContainsArea));
            EXPECT_FALSE(intersects(polygon, area, ContainedByArea));

            // only intersects
            area = S2Cap(S2LatLng::FromDegrees(1.1, 1.1).Normalized().ToPoint(), S1ChordAngle::Degrees(0.15));
            EXPECT_FALSE(intersects(polygon, area, ContainsCenter));
            EXPECT_TRUE(intersects(polygon, area, IntersectsArea));
            EXPECT_FALSE(intersects(polygon, area, ContainsArea));
            EXPECT_FALSE(intersects(polygon, area, ContainedByArea));

            // polygon contains area
            area = S2Cap(S2LatLng::FromDegrees(0.6, 0.6).Normalized().ToPoint(), S1ChordAngle::Degrees(0.15));
            EXPECT_TRUE(intersects(polygon, area, ContainsCenter));
            EXPECT_TRUE(intersects(polygon, area, IntersectsArea));
            EXPECT_TRUE(intersects(polygon, area, ContainsArea));
            EXPECT_FALSE(intersects(polygon, area, ContainedByArea));

            // polygon contained by area
            area = S2Cap(S2LatLng::FromDegrees(0.6, 0.6).Normalized().ToPoint(), S1ChordAngle::Degrees(0.9));
            EXPECT_TRUE(intersects(polygon, area, ContainsCenter));
            EXPECT_TRUE(intersects(polygon, area, IntersectsArea));
            EXPECT_FALSE(intersects(polygon, area, ContainsArea));
            EXPECT_TRUE(intersects(polygon, area, ContainedByArea));
        }

        TEST(SBSearchTests, PolygonIntersectsArea)
        {
            S2Polygon polygon, area;

            // do not intersect
            util::make_polygon(util::make_vertices("0:0, 0:1, 1:1, 1:0"), polygon);
            util::make_polygon(util::make_vertices("1.1:1.1, 1.1:2.1, 2.1:2.1, 2.1:1.1"), area);
            EXPECT_FALSE(intersects(polygon, area, ContainsCenter));
            EXPECT_FALSE(intersects(polygon, area, IntersectsArea));
            EXPECT_FALSE(intersects(polygon, area, ContainsArea));
            EXPECT_FALSE(intersects(polygon, area, ContainedByArea));

            // only intersects
            util::make_polygon(util::make_vertices("0.9:0.9, 0.9:2.1, 2.1:2.1, 2.1:0.9"), area);
            EXPECT_FALSE(intersects(polygon, area, ContainsCenter));
            EXPECT_TRUE(intersects(polygon, area, IntersectsArea));
            EXPECT_FALSE(intersects(polygon, area, ContainsArea));
            EXPECT_FALSE(intersects(polygon, area, ContainedByArea));

            // polygon contains area
            util::make_polygon(util::make_vertices("0.9:0.9, 0.9:0.95, 0.95:0.95, 0.95:0.9"), area);
            EXPECT_TRUE(intersects(polygon, area, ContainsCenter));
            EXPECT_TRUE(intersects(polygon, area, IntersectsArea));
            EXPECT_TRUE(intersects(polygon, area, ContainsArea));
            EXPECT_FALSE(intersects(polygon, area, ContainedByArea));

            // polygon contained by area
            util::make_polygon(util::make_vertices("-1:-1, -1:2, 2:2, 2:-1"), area);
            EXPECT_TRUE(intersects(polygon, area, ContainsCenter));
            EXPECT_TRUE(intersects(polygon, area, IntersectsArea));
            EXPECT_FALSE(intersects(polygon, area, ContainsArea));
            EXPECT_TRUE(intersects(polygon, area, ContainedByArea));
        }
    }
}