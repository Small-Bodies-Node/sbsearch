#include <gtest/gtest.h>

#include "../ephemeris.h"
#include "../exceptions.h"
#include "../moving_target.h"
#include "add.h"
#include "get.h"
#include "remove.h"
#include "sbsdb_test.h"

using namespace sbsearch;
using namespace sbsearch::sbsdb;

namespace sbsearch::sbsdb::testing
{
    // SBSearchDatabaseTest fixture defined in sbsdb_test.cc

    TEST_F(SBSearchDatabaseTest, EphemerisIO)
    {
        MovingTarget encke{"2P"};
        Ephemeris eph{encke,
                      {{0, 10, 1, 0, 1, 0.1, 90, 0, 1, 180, 0, 0, 0, 10, -1},
                       {1, 11, 2, 0, 5, 0.5, 90, 1, 0, 0, 180, 30, 0, 20, 5},
                       {2, 12, 3, 0, 10, 1.0, 90, 2, 1, 90, 80, 90, 0, 30, 10}}};

        // The target is not in the database, so we expect an error
        EXPECT_THROW(add::ephemeris(&db, eph), MovingTargetError);

        // Add the target, verify that the id was updated
        add::moving_target(&db, encke);
        EXPECT_NE(encke.moving_target_id(), eph.target().moving_target_id());

        // Fix the target, and then we can add the ephemeris data
        eph.target(encke);
        add::ephemeris(&db, eph);

        // Get the data back
        Ephemeris test;
        test = get::ephemeris(&db, eph.target());
        EXPECT_EQ(test, eph);

        // Get a subset of data
        test = get::ephemeris(&db, eph.target(), 0.5, 1.5);
        EXPECT_EQ(test, eph[1]);

        // This target does not match database copy
        MovingTarget wrong_id{"1P", eph.target().moving_target_id()};
        Ephemeris other{wrong_id, eph.data()};
        EXPECT_THROW(add::ephemeris(&db, other), MovingTargetError);

        // Remove some data
        remove::ephemeris(&db, eph.target(), 1.5, 10);
        test = get::ephemeris(&db, eph.target());
        EXPECT_NE(test, eph);
        EXPECT_EQ(test, eph.slice(0, 2));

        // Remove all
        remove::ephemeris(&db, eph.target());
        test = get::ephemeris(&db, eph.target());
        EXPECT_EQ(test.num_vertices(), 0);
    }

    TEST_F(SBSearchDatabaseTest, EphemerisDateRange)
    {
        MovingTarget encke{"2P"};
        EXPECT_THROW(get::ephemeris_date_range(&db, encke), MovingTargetError); // need obsid

        add::moving_target(&db, encke);
        auto range = get::ephemeris_date_range(&db, encke);
        EXPECT_FALSE(range.first); // no ephemeris data
        EXPECT_FALSE(range.second);

        Ephemeris encke_eph{encke,
                            {{0, 10, 1, 0, 1, 0.1, 90, 0, 1, 180, 0, 0, 0, 10, -1},
                             {1, 11, 2, 0, 5, 0.5, 90, 1, 0, 0, 180, 30, 0, 20, 5},
                             {2, 12, 3, 0, 10, 1.0, 90, 2, 1, 90, 80, 90, 0, 30, 10}}};
        add::ephemeris(&db, encke_eph);

        MovingTarget tempel1("9P");
        add::moving_target(&db, tempel1);
        Ephemeris tempel1_eph{tempel1,
                              {{1.5, 10, 1, 0, 1, 0.1, 90, 0, 1, 180, 0, 0, 0, 10, -1},
                               {2.5, 11, 2, 0, 5, 0.5, 90, 1, 0, 0, 180, 30, 0, 20, 5},
                               {3.5, 12, 3, 0, 10, 1.0, 90, 2, 1, 90, 80, 90, 0, 30, 10}}};
        add::ephemeris(&db, tempel1_eph);

        range = get::ephemeris_date_range(&db, encke);
        EXPECT_EQ(range.first, 0);
        EXPECT_EQ(range.second, 2);

        range = get::ephemeris_date_range(&db, tempel1);
        EXPECT_EQ(range.first, 1.5);
        EXPECT_EQ(range.second, 3.5);
    }
}