#include "config.h"

#include <cmath>
#include <vector>
#include <gtest/gtest.h>

#include "math.h"

namespace sbsearch::util
{
    TEST(UtilMathTests, IsIncreasing)
    {
        EXPECT_TRUE(is_increasing({1, 2, 3, 5, 100}));
        EXPECT_FALSE(is_increasing({1, 20, 3, 5, 100}));
        EXPECT_FALSE(is_increasing({1, 2, 3, 3, 100}));
        EXPECT_TRUE(is_increasing({0}));
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
