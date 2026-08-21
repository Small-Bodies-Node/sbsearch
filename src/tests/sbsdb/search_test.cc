#include <unordered_set>
#include <optional>
#include <gtest/gtest.h>

#include "observation.h"
#include "sbsdb/add.h"
#include "sbsdb/search.h"
#include "./sbsdb_test.h"

using namespace sbsearch;
using namespace sbsearch::sbsdb;

namespace sbsearch::sbsdb::testing
{
    // SBSearchDatabaseTest fixture defined in sbsdb_test.cc

    TEST_F(SBSearchDatabaseTest, SearchObservations)
    {
        Observations obs({{"test source", "X05", "a", 0, 1, "0:0, 0:1, 1:1", {"a", "b", "c"}, {}, "b"},
                          {"test source", "X05", "b", 1, 2, "0:0, 0:1, 1:1", {"b", "c", "d"}, {}, "c"},
                          {"test source", "X05", "c", 2, 3, "0:0, 0:1, 1:1", {"c", "d", "e"}, {}, "d"},
                          {"another test source", "T05", "d", 4, 5, "0:0, 0:1, 1:1", {"d", "e", "f"}, {}, "e"}});
        add::observations(&db, obs);

        // search for observations matching term a
        int n = search::observations(&db, {"a"});
        Observations matches = search::results(&db);
        EXPECT_EQ(n, 1);
        EXPECT_EQ(matches.size(), 1);

        // a or f
        n = search::observations(&db, {"a", "f"});
        matches = search::results(&db);
        EXPECT_EQ(n, 2);
        EXPECT_EQ(matches.size(), 2);

        // c or f
        n = search::observations(&db, {"c", "f"});
        matches = search::results(&db);
        EXPECT_EQ(n, 4);
        EXPECT_EQ(matches.size(), 4);

        // g
        n = search::observations(&db, {"g"});
        matches = search::results(&db);
        EXPECT_EQ(n, 0);
        EXPECT_EQ(matches.size(), 0);

        // test observation time limits
        // start
        n = search::observations(&db, {"e"}, {.mjd_start = 2});
        matches = search::results(&db);
        EXPECT_EQ(n, 2);
        EXPECT_EQ(matches.size(), 2);

        n = search::observations(&db, {"e"}, {.mjd_start = 3.5});
        matches = search::results(&db);
        EXPECT_EQ(n, 1);
        EXPECT_EQ(matches.size(), 1);

        // stop
        n = search::observations(&db, {"e"}, {.mjd_stop = 1});
        matches = search::results(&db);
        EXPECT_EQ(n, 0);
        EXPECT_EQ(matches.size(), 0);

        n = search::observations(&db, {"e"}, {.mjd_stop = 3});
        matches = search::results(&db);
        EXPECT_EQ(n, 1);
        EXPECT_EQ(matches.size(), 1);

        n = search::observations(&db, {"e"}, {.mjd_stop = 5});
        matches = search::results(&db);
        EXPECT_EQ(n, 2);
        EXPECT_EQ(matches.size(), 2);

        // start-stop
        n = search::observations(&db, {"e"}, {.mjd_start = 2, .mjd_stop = 2.5});
        matches = search::results(&db);
        EXPECT_EQ(n, 0);
        EXPECT_EQ(matches.size(), 0);

        n = search::observations(&db, {"e"}, {.mjd_start = 2, .mjd_stop = 3});
        matches = search::results(&db);
        EXPECT_EQ(n, 1);
        EXPECT_EQ(matches.size(), 1);

        n = search::observations(&db, {"e"}, {.mjd_start = 2.5, .mjd_stop = 4.5});
        matches = search::results(&db);
        EXPECT_EQ(n, 0);
        EXPECT_EQ(matches.size(), 0);

        n = search::observations(&db, {"e"}, {.mjd_start = 3, .mjd_stop = 5});
        matches = search::results(&db);
        EXPECT_EQ(n, 1);
        EXPECT_EQ(matches.size(), 1);

        // search by source
        n = search::observations(&db, {"b", "e"}, {.source = "test source"});
        matches = search::results(&db);
        EXPECT_EQ(n, 3);
        EXPECT_EQ(matches.size(), 3);

        n = search::observations(&db, {"b", "e"}, {.source = "another test source"});
        matches = search::results(&db);
        EXPECT_EQ(n, 1);
        EXPECT_EQ(matches.size(), 1);

        // test appending multiple searches
        n = search::observations(&db, {"a"});
        EXPECT_EQ(n, 1);

        n = search::observations(&db, {"c"});
        EXPECT_EQ(n, 3);

        n = search::observations(&db, {"f"});
        matches = search::results(&db);
        EXPECT_EQ(n, 4);
        EXPECT_EQ(matches.size(), 4);

        // temporary results table does not exist
        EXPECT_THROW(search::results(&db), std::exception);
    }
}