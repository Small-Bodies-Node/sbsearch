#include "config.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <vector>
#include <gtest/gtest.h>

#include "util/math.h"

using std::vector;

namespace sbsearch::util
{
    TEST(UtilMathTests, IsMonotonicallyIncreasing)
    {
        EXPECT_TRUE(is_monotonically_increasing(vector<double>({1, 2, 3, 5, 100})));
        EXPECT_TRUE(is_monotonically_increasing(vector<double>({1, 2, 3, 3, 100})));
        EXPECT_TRUE(is_monotonically_increasing(vector<double>({0})));

        EXPECT_FALSE(is_monotonically_increasing(vector<double>({1, 20, 3, 5, 100})));
        EXPECT_FALSE(is_monotonically_increasing(vector<std::optional<double>>({1, 2, 3, std::nullopt, 100})));
    }

    TEST(UtilMathTests, IsStrictlyIncreasing)
    {
        EXPECT_TRUE(is_strictly_increasing(vector<double>({1, 2, 3, 5, 100})));
        EXPECT_TRUE(is_strictly_increasing(vector<double>({0})));

        EXPECT_FALSE(is_strictly_increasing(vector<double>({1, 2, 3, 3, 100})));
        EXPECT_FALSE(is_strictly_increasing(vector<double>({1, 20, 3, 5, 100})));
        EXPECT_FALSE(is_strictly_increasing(vector<std::optional<double>>({1, 2, 3, std::nullopt, 100})));
    }

    TEST(UtilMathTests, Interp)
    {
        EXPECT_EQ(interp(1, 2, 0.5), 1.5);
        EXPECT_EQ(interp(1, 2, 0.2), 1.2);
        EXPECT_EQ(interp(2, 1, 0.2), 1.8);
        EXPECT_EQ(interp(1, 2, 2), 3);
        EXPECT_EQ(interp(1, 2, 0), 1);
        EXPECT_EQ(interp(1, 2, -0.5), 0.5);
    }
}
