#include "config.h"

#include <set>
#include <string>
#include <vector>
#include <gtest/gtest.h>

#include "moving_target.h"

using std::set;
using std::string;
using std::vector;

namespace sbsearch
{
    namespace testing
    {
        TEST(MovingTargetTests, MovingTargetInit)
        {
            MovingTarget target;
            EXPECT_EQ(target.designation(), "");
            EXPECT_EQ(target.moving_target_id(), UNDEF_MOVING_TARGET_ID);
            EXPECT_EQ(target.alternate_names(), set<string>());
            EXPECT_EQ(target.small_body(), true);

            target = MovingTarget("2P");
            EXPECT_EQ(target.designation(), "2P");
            EXPECT_EQ(target.moving_target_id(), UNDEF_MOVING_TARGET_ID);
            EXPECT_EQ(target.alternate_names(), set<string>());
            EXPECT_EQ(target.small_body(), true);

            target = MovingTarget("2P", 1);
            EXPECT_EQ(target.designation(), "2P");
            EXPECT_EQ(target.moving_target_id(), 1);
            EXPECT_EQ(target.alternate_names(), set<string>());
            EXPECT_EQ(target.small_body(), true);

            MovingTarget new_target(target);
            EXPECT_EQ(target, new_target);
        }

        TEST(MovingTargetTests, MovingTargetEquality)
        {
            MovingTarget a("2P", 1);
            MovingTarget b("2P", 1);
            EXPECT_EQ(a, b);

            a.designation("1P");
            EXPECT_NE(a, b);
            a.designation("2P");
            EXPECT_EQ(a, b);

            a.moving_target_id(2);
            EXPECT_NE(a, b);
            a.moving_target_id(1);
            EXPECT_EQ(a, b);

            a.add_name("Encke");
            EXPECT_NE(a, b);
            b.add_name("Encke");
            EXPECT_EQ(a, b);
        }

        TEST(MovingTargetTests, MovingTargetDesignation)
        {
            MovingTarget a("P/2003 CC22", 1);

            // update designation, discarding the old one
            a.designation("452P");
            EXPECT_EQ(a.designation(), "452P");
            EXPECT_EQ(a.alternate_names(), set<string>());
        }

        TEST(MovingTargetTests, MovingTargetAddName)
        {
            MovingTarget a("P/2003 CC22", 1);
            EXPECT_EQ(a.alternate_names(), set<string>());

            a.add_name("Sheppard-Jewitt");
            EXPECT_EQ(a.alternate_names(), set<string>({"Sheppard-Jewitt"}));

            // update designation, keep old designation as an alternate name
            a.add_name("452P", true);
            EXPECT_EQ(a.alternate_names(), set<string>({"P/2003 CC22", "Sheppard-Jewitt"}));
        }

        TEST(MovingTargetTests, MovingTargetAddNames)
        {
            MovingTarget a("2P", 1);
            EXPECT_EQ(a.alternate_names(), set<string>());

            vector<string> names{"Encke", "2P/Encke"};
            a.add_names(names.begin(), names.end());
            EXPECT_EQ(a.alternate_names(), set<string>({"Encke", "2P/Encke"}));
        }

        TEST(MovingTargetTests, MovingTargetToString)
        {
            MovingTarget target("452P", 1);
            EXPECT_EQ(to_string(target), "452P (ID=1; small body=true)");

            target.add_name("Sheppard-Jewitt");
            target.add_name("P/2003 CC22");
            EXPECT_EQ(to_string(target), "452P (ID=1; P/2003 CC22, Sheppard-Jewitt; small body=true)");
        }
    }
}