#ifndef SBS_TESTS_ENV_H_
#define SBS_TESTS_ENV_H_

#include <string>
#include <gtest/gtest.h>

class Environment : public ::testing::Environment
{
public:
    Environment(int argc, char **argv)
    {
        for (int i = 0; i < argc; i++)
            if (std::string(argv[i]) == "--remote-data")
                skip_remote_ = false;
    }

    ~Environment() override {}

    // Override this to define how to set up the environment.
    void SetUp() override {}

    // Override this to define how to tear down the environment.
    void TearDown() override
    {
        if (skip_remote_)
            std::cout << "Remote tests were skipped.  Run them with --remote-data.\n";
    }

    inline const bool skip_remote()
    {
        return skip_remote_;
    }

private:
    bool skip_remote_ = true;
};

// Pointer to safely access the global environment in main()
extern Environment *env;

#endif
