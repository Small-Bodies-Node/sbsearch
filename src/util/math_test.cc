#include "config.h"

#include <cmath>
#include <optional>
#include <vector>
#include <gtest/gtest.h>

#include "math.h"

using std::vector;

namespace sbsearch::util
{
    TEST(UtilMathTests, IsIncreasing)
    {
        EXPECT_TRUE(is_increasing(vector<double>({1, 2, 3, 5, 100})));
        EXPECT_FALSE(is_increasing(vector<double>({1, 20, 3, 5, 100})));
        EXPECT_FALSE(is_increasing(vector<double>({1, 2, 3, 3, 100})));
        EXPECT_TRUE(is_increasing(vector<double>({0})));
        EXPECT_FALSE(is_increasing(vector<std::optional<double>>({1, 2, 3, std::nullopt, 100})));
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
