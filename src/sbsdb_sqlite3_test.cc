#include "sbsdb.h"
#include "sbsdb_sqlite3.h"
#include "observation.h"
#include "sbsearch.h"

#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string>

#include <gtest/gtest.h>
#include <s2/s2latlng.h>
#include <s2/s2latlng_rect.h>
#include <s2/s2metrics.h>
#include <s2/s2point.h>
#include <s2/s2region_term_indexer.h>
#include <sqlite3.h>

using sbsearch::Observation;
using sbsearch::SBSearchDatabase;
using sbsearch::SBSearchDatabaseSqlite3;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

class SBSSearchDatabaseSqlite3Test : public ::testing::Test
{
public:
    SBSearchDatabaseSqlite3 *sbsdb;
    SBSSearchDatabaseSqlite3Test()
    {
        SBSearchDatabaseSqlite3::Options options;
        options.max_spatial_cells(8);
        options.max_spatial_resolution(584.4);
        options.min_spatial_resolution(34.4);
        sbsdb = new SBSearchDatabaseSqlite3(":memory:", options);
        sbsdb->setup_tables();

        sbsdb->add_observation(Observation(59252.1, 59252.2, "1:3, 2:3, 2:4, 1:4"));
        sbsdb->add_observation(Observation(59252.21, 59252.31, "2:3, 3:3, 3:4, 1:4"));
    }
};

namespace sbsearch
{
    namespace testing
    {
        TEST(SBSearchDatabaseSqlite3Tests, SBSearchDatabaseSqlite3Init)
        {
            SBSearchDatabaseSqlite3 sbsdb(""); // open a temporary private file
            EXPECT_EQ(sbsdb.options().max_spatial_cells(), 8);
            EXPECT_EQ(sbsdb.options().max_spatial_level(), 12);
            EXPECT_EQ(sbsdb.options().min_spatial_level(), 4);
            sbsdb.close(); // close it
            EXPECT_THROW(sbsdb.drop_time_indices(), std::runtime_error);

            // try to open the root directory as a database file
            EXPECT_THROW(SBSearchDatabaseSqlite3("/"), std::runtime_error);

            // verify that options are correctly set
            SBSearchDatabase::Options options;
            options.max_spatial_cells(12);
            options.max_spatial_resolution(10);
            options.min_spatial_resolution(0.1);
            EXPECT_EQ(options.max_spatial_cells(), 12);
            EXPECT_EQ(options.max_spatial_level(), 16);
            EXPECT_EQ(options.min_spatial_level(), 10);
            SBSearchDatabaseSqlite3 sbsdb2(":memory:", options);
            EXPECT_EQ(sbsdb2.options().max_spatial_cells(), 12);
            EXPECT_EQ(options.max_spatial_level(), 16);
            EXPECT_EQ(options.min_spatial_level(), 10);
        }

        TEST(SBSearchDatabaseSqlite3Tests, SBSearchDatabaseSqlite3SetupTables)
        {
            SBSearchDatabaseSqlite3 sbsdb(":memory:");
            Observation obs(0, 1, "0:0, 0:1, 1:1");
            EXPECT_THROW(sbsdb.add_observation(obs), std::runtime_error);
            sbsdb.setup_tables();
            EXPECT_NO_THROW(sbsdb.add_observation(obs));
        }

        TEST(SBSearchDatabaseSqlite3Tests, SBSearchDatabaseSqlite3GetOneValue)
        {
            SBSearchDatabaseSqlite3 sbsdb(":memory:");
            sbsdb.setup_tables();
            Observation obs(0, 1, "0:0, 0:1, 1:1");
            sbsdb.add_observation(obs);
            double mjd = sbsdb.get_one_value<double>("SELECT mjd_start FROM observations LIMIT 1");
            EXPECT_EQ(mjd, 0);
        }
    }
}
