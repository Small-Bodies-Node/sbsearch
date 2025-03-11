#include "config.h"

#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <gtest/gtest.h>

#include "ephemeris.h"
#include "exceptions.h"
#include "indexer.h"
#include "moving_target.h"
#include "observation.h"
#include "observatory.h"
#include "sbsdb.h"
#include "sbsdb_sqlite3.h"

using namespace sbsearch;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace testing
{
    TEST(SBSearchDatabaseSqlite3Tests, Init)
    {
        // open a temporary private file
        SBSearchDatabaseSqlite3 sbsdb("");
        sbsdb.close();

        // try to open the root directory as a database file
        EXPECT_THROW(SBSearchDatabaseSqlite3("/"), std::runtime_error);
    }

    TEST(SBSearchDatabaseSqlite3Tests, DropCreateObservationsIndices)
    {
        SBSearchDatabaseSqlite3 sbsdb(":memory:");
        sbsdb.setup_tables();
        EXPECT_EQ(sbsdb.get_int("SELECT COUNT(*) FROM sqlite_master WHERE type='index' and name='idx_observations_mjd_start';"), 1);
        EXPECT_EQ(sbsdb.get_int("SELECT COUNT(*) FROM sqlite_master WHERE type='index' and name='idx_observations_mjd_stop';"), 1);
        sbsdb.drop_observations_indices();
        EXPECT_EQ(sbsdb.get_int("SELECT COUNT(*) FROM sqlite_master WHERE type='index' and name='idx_observations_mjd_start';"), 0);
        EXPECT_EQ(sbsdb.get_int("SELECT COUNT(*) FROM sqlite_master WHERE type='index' and name='idx_observations_mjd_stop';"), 0);
        sbsdb.create_observations_indices();
        EXPECT_EQ(sbsdb.get_int("SELECT COUNT(*) FROM sqlite_master WHERE type='index' and name='idx_observations_mjd_start';"), 1);
        EXPECT_EQ(sbsdb.get_int("SELECT COUNT(*) FROM sqlite_master WHERE type='index' and name='idx_observations_mjd_stop';"), 1);
    }

    TEST(SBSearchDatabaseSqlite3Tests, ErrorIfClosed)
    {
        SBSearchDatabaseSqlite3 sbsdb(":memory:");
        sbsdb.close();
        EXPECT_THROW(sbsdb.execute_sql("SELECT 1"), std::runtime_error);
    }
}