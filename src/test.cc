#include "gtest/gtest.h"

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    // remote tests must be explicitly specified with --gtest_filter=*:*Remote*
    if (GTEST_FLAG_GET(filter) == "*")
        GTEST_FLAG_SET(filter, "-*Remote*");
    return RUN_ALL_TESTS();
}
