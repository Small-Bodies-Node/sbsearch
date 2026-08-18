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

#include "polygons.h"
#include "util/string.h"

using std::ceil;
using std::floor;
using std::string;
using std::vector;

namespace sbsearch::util
{
    TEST(UtilPolygonTests, MakePolygonSimple)
    {
        S2Polygon polygon;
        make_polygon_simple(make_vertices("0:0, 1:0, 1:1, 0:1"), polygon);
        // note: s2geometry's text format is lat:lng
        std::unique_ptr<S2Polygon> expected = s2textformat::MakePolygonOrDie("0:0, 0:1, 1:1, 1:0");
        EXPECT_TRUE(polygon.BoundaryEquals(*expected.get()));
    }

    TEST(UtilPolygonTests, MakePolygon)
    {
        S2Polygon polygon;
        // order does not matter
        make_polygon(make_vertices("0:1, 1:1, 1:0, 0:0"), polygon);
        // note: s2geometry's text format is lat:lng
        std::unique_ptr<S2Polygon> expected = s2textformat::MakePolygonOrDie("0:0, 0:1, 1:1, 1:0");
        EXPECT_TRUE(polygon.BoundaryEquals(*expected.get()));
    }

    TEST(UtilPolygonTests, PaddedPolygon)
    {
        // results based on diagonstics.cc

        S2Polygon polygon, result;
        make_polygon(make_vertices("0:0, 1:0, 1:1"), polygon);
        padded_polygon(polygon, 0.2 * 60, result);
        EXPECT_EQ(format_vertices(result), "0.000000:-0.200000, 1.000000:-0.200000, 1.078008:-0.184160, 1.143658:-0.139148, 1.186553:-0.072096, 1.200000:0.000000, 1.200030:0.999994, 1.184192:1.078002, 1.139176:1.143655, 1.072111:1.186553, 0.993622:1.199898, 0.916143:1.181578, 0.858561:1.141429, -0.141432:0.141410, -0.185386:0.075047, -0.199974:-0.003204, -0.182887:-0.080947, -0.136830:-0.145869, -0.069098:-0.187684");

        padded_polygon(polygon, 2.0 * 60, result);
        EXPECT_EQ(format_vertices(result), "0.000000:-2.000000, 1.000000:-2.000000, 1.780631:-1.841419, 2.437301:-1.390889, 2.865971:-0.719955, 3.000000:0.000000, 3.000304:0.999391, 2.842308:1.779711, 2.392148:2.436582, 1.720856:2.865745, 0.934903:2.998942, 0.159309:2.814960, -0.415218:2.413873, -1.414608:1.413962, -1.854083:0.750049, -1.999733:-0.032661, -1.828606:-0.810198, -1.367674:-1.459408, -0.689756:-1.877340");
    }
}
