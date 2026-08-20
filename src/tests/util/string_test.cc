#include "config.h"

#include <stdexcept>
#include <string>

#include <gtest/gtest.h>
#include <s2/s2latlng.h>
#include <s2/s2latlng_rect.h>
#include <s2/s2point.h>

#include "util/string.h"

using std::string;
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
}
