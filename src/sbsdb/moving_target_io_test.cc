#include <string>
#include <vector>
#include <gtest/gtest.h>

#include "../exceptions.h"
#include "../moving_target.h"
#include "../observation.h"
#include "../observatory.h"
#include "add.h"
#include "get.h"
#include "remove.h"
#include "sbsdb_test.h"
#include "update.h"

using namespace sbsearch;
using namespace sbsearch::sbsdb;
using std::string;
using std::vector;

namespace sbsearch::sbsdb::testing
{
    // SBSearchDatabaseTest fixture defined in sbsdb_test.cc

    TEST_F(SBSearchDatabaseTest, MovingTargetIO)
    {
        MovingTarget encke("2P");
        MovingTarget ceres("1");
        MovingTarget mercury("1", false);

        // verify encke should fail: not yet in the database
        EXPECT_THROW(verify::moving_target(&db, encke), MovingTargetError);

        // add to the database, expect an updated moving_target_id
        EXPECT_FALSE(encke.moving_target_id());
        add::moving_target(&db, encke);
        EXPECT_EQ(encke.moving_target_id(), 1);

        // verify encke should now pass
        EXPECT_NO_THROW(verify::moving_target(&db, encke));

        // add another, it should be 2
        add::moving_target(&db, ceres);
        EXPECT_EQ(ceres.moving_target_id(), 2);

        // and this should be 3
        add::moving_target(&db, mercury);
        EXPECT_EQ(mercury.moving_target_id(), 3);

        // get them from the database
        MovingTarget test;
        test = get::moving_target(&db, 1);
        EXPECT_EQ(test.designation(), encke.designation());
        EXPECT_EQ(test.moving_target_id(), encke.moving_target_id());
        EXPECT_EQ(test.alternate_names(), encke.alternate_names());
        EXPECT_EQ(test.small_body(), encke.small_body());

        test = get::moving_target(&db, "2P");
        EXPECT_EQ(test.designation(), encke.designation());
        EXPECT_EQ(test.moving_target_id(), encke.moving_target_id());
        EXPECT_EQ(test.alternate_names(), encke.alternate_names());
        EXPECT_EQ(test.small_body(), encke.small_body());

        test = get::moving_target(&db, "1");
        EXPECT_EQ(test.designation(), ceres.designation());
        EXPECT_EQ(test.moving_target_id(), ceres.moving_target_id());
        EXPECT_EQ(test.alternate_names(), ceres.alternate_names());
        EXPECT_EQ(test.small_body(), ceres.small_body());

        test = get::moving_target(&db, 2);
        EXPECT_EQ(test.designation(), ceres.designation());
        EXPECT_EQ(test.moving_target_id(), ceres.moving_target_id());
        EXPECT_EQ(test.alternate_names(), ceres.alternate_names());
        EXPECT_EQ(test.small_body(), ceres.small_body());

        test = get::moving_target(&db, "1", false);
        EXPECT_EQ(test.designation(), mercury.designation());
        EXPECT_EQ(test.moving_target_id(), mercury.moving_target_id());
        EXPECT_EQ(test.alternate_names(), mercury.alternate_names());
        EXPECT_EQ(test.small_body(), mercury.small_body());

        test = get::moving_target(&db, 3);
        EXPECT_EQ(test.designation(), mercury.designation());
        EXPECT_EQ(test.moving_target_id(), mercury.moving_target_id());
        EXPECT_EQ(test.alternate_names(), mercury.alternate_names());
        EXPECT_EQ(test.small_body(), mercury.small_body());

        // add an alternate name to encke and it no longer matches the db copy
        encke.add_name("2P/Encke");
        EXPECT_THROW(verify::moving_target(&db, encke), MovingTargetError);

        // update and it should pass
        update::moving_target(&db, encke);
        EXPECT_NO_THROW(verify::moving_target(&db, encke));

        // try getting Encke via alt name
        test = get::moving_target(&db, "2P/Encke");
        EXPECT_EQ(test.designation(), encke.designation());
        EXPECT_EQ(test.moving_target_id(), encke.moving_target_id());
        EXPECT_EQ(test.alternate_names(), encke.alternate_names());

        // add a few names to Ceres
        vector<string> names{"(1) Ceres", "Ceres", "A801 AA"};
        ceres.add_names(names.begin(), names.end());
        update::moving_target(&db, ceres);
        test = get::moving_target(&db, "A801 AA");
        EXPECT_EQ(test.designation(), ceres.designation());
        EXPECT_EQ(test.moving_target_id(), ceres.moving_target_id());
        EXPECT_EQ(test.alternate_names(), ceres.alternate_names());
        EXPECT_EQ(test.alternate_names().size(), 3);

        // and Mercury
        names = vector<string>{"Hermes", "Nabu"};
        mercury.add_names(names.begin(), names.end());
        update::moving_target(&db, mercury);
        test = get::moving_target(&db, "Nabu", false);
        EXPECT_EQ(test.designation(), mercury.designation());
        EXPECT_EQ(test.moving_target_id(), mercury.moving_target_id());
        EXPECT_EQ(test.alternate_names(), mercury.alternate_names());
        EXPECT_EQ(test.alternate_names().size(), 2);

        // Try to add an object that already exists.
        MovingTarget halley("1P");
        halley.moving_target_id(1);
        EXPECT_THROW(add::moving_target(&db, halley), MovingTargetError);

        MovingTarget duplicate_ceres("1");
        EXPECT_THROW(add::moving_target(&db, duplicate_ceres), MovingTargetError);

        // Removing an object that does not exist is just a warning.
        MovingTarget new_comet("1000P");
        new_comet.moving_target_id(9123);
        remove::moving_target(&db, new_comet);

        // Get objects that do not exist
        EXPECT_THROW(get::moving_target(&db, 1000), MovingTargetError);
        EXPECT_EQ(get::moving_target(&db, "asdf"), MovingTarget("asdf"));
        EXPECT_FALSE(get::moving_target(&db, "Nabu", true).moving_target_id());
        EXPECT_FALSE(get::moving_target(&db, "Ceres", false).moving_target_id());
    }

}