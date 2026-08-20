#include <string>
#include <gtest/gtest.h>
#include "tests/env.h"

// Define global environment pointer
Environment *env = nullptr;

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    env = new Environment(argc, argv);
    ::testing::AddGlobalTestEnvironment(env);

    return RUN_ALL_TESTS();
}
