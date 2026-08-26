#include <stdexcept>
#include <string>
#include <string_view>
#include <gtest/gtest.h>
#include <s2/s2latlng.h>
#include <s2/s2latlng_rect.h>
#include <s2/s2point.h>

#include "config.h"
#include "util/string.h"

using std::string;
using std::string_view;
using std::vector;
using namespace std::literals::string_view_literals;

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

    TEST(UtilStringTests, Strip)
    {
        EXPECT_EQ(strip(""), ""sv);
        EXPECT_EQ(strip(" "), ""sv);
        EXPECT_EQ(strip("   "), ""sv);
        EXPECT_EQ(strip("asdf"), "asdf"sv);
        EXPECT_EQ(strip(" asdf "), "asdf"sv);
        EXPECT_EQ(strip("  asdf  "), "asdf"sv);
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

    TEST(UtilStringTests, StringViewToDouble)
    {
        EXPECT_FLOAT_EQ(svtod("1.234"sv), 1.234);
        EXPECT_FLOAT_EQ(svtod("1.2345"sv), 1.2345);
        EXPECT_FLOAT_EQ(svtod("1.23456"sv), 1.23456);
        EXPECT_FLOAT_EQ(svtod("1.234567"sv), 1.234567);
        EXPECT_FLOAT_EQ(svtod("1.2345678"sv), 1.2345678);

        EXPECT_THROW(svtod("asdf"sv), std::invalid_argument);
    }

    TEST(UtilStringTests, DoubleToString)
    {
        EXPECT_EQ(dtos(1.2345678, 3), "1.235");
        EXPECT_EQ(dtos(1.2345678, 4), "1.2346");
        EXPECT_EQ(dtos(1.2345678, 5), "1.23457");
        EXPECT_EQ(dtos(1.2345678, 6), "1.234568");
        EXPECT_EQ(dtos(1.2345678, 7), "1.2345678");
        EXPECT_EQ(dtos(1.2345678), "1.2345678");
        EXPECT_EQ(dtos(1.2345678, 8), "1.23456780");
    }
}
