#include <unordered_set>
#include <gtest/gtest.h>

#include "../observation.h"
#include "add.h"
#include "find.h"
#include "sbsdb_test.h"

using namespace sbsearch;
using namespace sbsearch::sbsdb;

namespace sbsearch::sbsdb::testing
{
    // SBSearchDatabaseTest fixture defined in sbsdb_test.cc

    TEST_F(SBSearchDatabaseTest, FindObservations)
    {
        Observations obs({{"test source", "X05", "a", 0, 1, "0:0, 0:1, 1:1", "a b c"},
                          {"test source", "X05", "b", 1, 2, "0:0, 0:1, 1:1", "b c d"},
                          {"test source", "X05", "c", 2, 3, "0:0, 0:1, 1:1", "c d e"},
                          {"another test source", "T05", "d", 4, 5, "0:0, 0:1, 1:1", "d e f"}});
        add::observations(&db, obs);

        // find observations matching term a
        std::unordered_set<Observation> matches;
        matches = find::observations(&db, {"a"});
        EXPECT_EQ(matches.size(), 1);

        // a or f
        matches = find::observations(&db, {"a", "f"});
        EXPECT_EQ(matches.size(), 2);

        // c or f
        matches = find::observations(&db, {"c", "f"});
        EXPECT_EQ(matches.size(), 4);

        // g
        matches = find::observations(&db, {"g"});
        EXPECT_EQ(matches.size(), 0);

        // test observation time limits
        // start
        matches = find::observations(&db, {"e"}, {.mjd_start = 2});
        EXPECT_EQ(matches.size(), 2);

        matches = find::observations(&db, {"e"}, {.mjd_start = 3.5});
        EXPECT_EQ(matches.size(), 1);

        // stop
        matches = find::observations(&db, {"e"}, {.mjd_stop = 1});
        EXPECT_EQ(matches.size(), 0);

        matches = find::observations(&db, {"e"}, {.mjd_stop = 3});
        EXPECT_EQ(matches.size(), 1);

        matches = find::observations(&db, {"e"}, {.mjd_stop = 5});
        EXPECT_EQ(matches.size(), 2);

        // start-stop
        matches = find::observations(&db, {"e"}, {.mjd_start = 2, .mjd_stop = 2.5});
        EXPECT_EQ(matches.size(), 0);

        matches = find::observations(&db, {"e"}, {.mjd_start = 2, .mjd_stop = 3});
        EXPECT_EQ(matches.size(), 1);

        matches = find::observations(&db, {"e"}, {.mjd_start = 2.5, .mjd_stop = 4.5});
        EXPECT_EQ(matches.size(), 0);

        matches = find::observations(&db, {"e"}, {.mjd_start = 3, .mjd_stop = 5});
        EXPECT_EQ(matches.size(), 1);

        // search by source
        matches = find::observations(&db, {"b", "e"}, {.source = "test source"});
        EXPECT_EQ(matches.size(), 3);

        matches = find::observations(&db, {"b", "e"}, {.source = "another test source"});
        EXPECT_EQ(matches.size(), 1);
    }

}