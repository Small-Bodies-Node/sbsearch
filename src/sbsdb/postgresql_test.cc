#include <optional>
#include <gtest/gtest.h>

#include "postgresql.h"
#include "sbsdb_test.h"

using namespace sbsearch;
using namespace sbsearch::sbsdb;
using std::string;

namespace sbsearch::sbsdb::testing
{
    // SBSearchDatabaseTest fixture defined in sbsdb_test.cc

    TEST_F(SBSearchDatabaseTest, SetupTables)
    {
        EXPECT_EQ(db.get_one<int>("SELECT COUNT(*) FROM configuration"), 4);
        EXPECT_EQ(db.get_one<int>("SELECT COUNT(*) FROM ephemerides"), 0);
        EXPECT_EQ(db.get_one<int>("SELECT COUNT(*) FROM found"), 0);
        EXPECT_EQ(db.get_one<int>("SELECT COUNT(*) FROM moving_targets"), 0);
        EXPECT_EQ(db.get_one<int>("SELECT COUNT(*) FROM observations"), 0);
        EXPECT_EQ(db.get_one<int>("SELECT COUNT(*) FROM observatories"), 0);

        EXPECT_NO_THROW(db.setup_tables());
    }

    TEST_F(SBSearchDatabaseTest, SelectFromMissingTable)
    {
        // try to get a value from a table that does not exist
        EXPECT_THROW(db.get_one<int>("SELECT observation_id FROM invalid_table LIMIT 1"),
                     std::runtime_error);
    }

    TEST_F(SBSearchDatabaseTest, GetInt)
    {
        auto value = db.get_one<int>("SELECT 1");
        EXPECT_EQ(value, 1);

        auto optional_value = db.get_one<optional<int>>("SELECT 1");
        EXPECT_EQ(optional_value, 1);

        optional_value = db.get_one<optional<int>>("SELECT NULL");
        EXPECT_FALSE(optional_value);
    }

    TEST_F(SBSearchDatabaseTest, GetInt64)
    {
        auto value = db.get_one<int64_t>("SELECT 1234567890123456789");
        EXPECT_EQ(value, 1234567890123456789);

        auto optional_value = db.get_one<optional<int64_t>>("SELECT 1234567890123456789");
        EXPECT_EQ(optional_value, 1234567890123456789);

        optional_value = db.get_one<optional<int64_t>>("SELECT NULL");
        EXPECT_FALSE(optional_value);
    }

    TEST_F(SBSearchDatabaseTest, GetDouble)
    {
        auto value = db.get_one<double>("SELECT 1.123");
        EXPECT_EQ(value, 1.123);

        auto optional_value = db.get_one<optional<double>>("SELECT 1.123");
        EXPECT_EQ(optional_value, 1.123);

        optional_value = db.get_one<optional<double>>("SELECT NULL");
        EXPECT_FALSE(optional_value);
    }

    TEST_F(SBSearchDatabaseTest, GetString)
    {
        auto value = db.get_one<string>("SELECT 'asdf'");
        EXPECT_EQ(value, "asdf");

        auto optional_value = db.get_one<optional<string>>("SELECT 'asdf'");
        EXPECT_EQ(optional_value, "asdf");

        optional_value = db.get_one<optional<string>>("SELECT NULL");
        EXPECT_FALSE(optional_value);
    }
}