#include <optional>
#include <gtest/gtest.h>

#include "../ephemeris.h"
#include "../moving_target.h"
#include "../observation.h"
#include "add.h"
#include "get.h"
#include "remove.h"
#include "sbsdb_test.h"

using namespace sbsearch;
using namespace sbsearch::sbsdb;

namespace sbsearch::sbsdb::testing
{
    // SBSearchDatabaseTest fixture defined in sbsdb_test.cc

    TEST_F(SBSearchDatabaseTest, FoundIO)
    {
        Observations observations({{"test source", "X05", "a", 0, 1, "0:0, 0:1, 1:1", "a b c", std::nullopt, "b"},
                                   {"test source", "X05", "b", 1, 2, "0:0, 0:1, 1:1", "b c d", std::nullopt, "c"}});
        add::observations(&db, observations);

        MovingTarget encke("2P");
        add::moving_target(&db, encke);

        Ephemeris eph{encke,
                      {{0, 10, 1, 0, 1, 0.1, 90, 0, 1, 180, 0, 0, 0, 10, -1},
                       {1, 11, 2, 0, 5, 0.5, 90, 1, 0, 0, 180, 30, 0, 20, 5},
                       {2, 12, 3, 0, 10, 1.0, 90, 2, 1, 90, 80, 90, 0, 30, 10}}};

        // these may not make sense, but it doesn't matter
        Founds founds;
        founds.append(Found(observations[0], eph.segment(0)));
        founds.append(Found(observations[1], eph.segment(1)));
        add::found(&db, founds);

        founds = get::found(&db, observations[0]);
        EXPECT_EQ(founds.size(), 1);
        EXPECT_EQ(founds[0], Found(observations[0], eph[0]));

        founds = get::found(&db, observations[1]);
        EXPECT_EQ(founds.size(), 1);
        EXPECT_EQ(founds[0], Found(observations[1], eph[1]));

        founds = get::found(&db, encke);
        EXPECT_EQ(founds.size(), 2);
        EXPECT_EQ(std::count(founds.begin(), founds.end(), Found(observations[0], eph[0])), 1);
        EXPECT_EQ(std::count(founds.begin(), founds.end(), Found(observations[1], eph[1])), 1);

        remove::found(&db, founds);
        founds = get::found(&db, observations[0]);
        EXPECT_EQ(founds.size(), 0);
        founds = get::found(&db, observations[1]);
        EXPECT_EQ(founds.size(), 0);
        founds = get::found(&db, encke);
        EXPECT_EQ(founds.size(), 0);
    }

}