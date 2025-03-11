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
#include "sbsdb_postgresql.h"

using namespace sbsearch;
using std::cerr;
using std::cout;
using std::string;
using std::vector;

class SBSearchDatabasePostgreSQLTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        sbsdb.execute_sql("BEGIN");
        sbsdb.execute_sql("SAVEPOINT test_setup");
    }

    void TearDown() override
    {
        sbsdb.execute_sql("ROLLBACK TO SAVEPOINT test_setup");
        sbsdb.execute_sql("COMMIT");
    }

    SBSearchDatabasePostgreSQL sbsdb{"dbname=sbsearch_test"};
};

namespace testing
{
    TEST(SBSearchDatabasePostgreSQLTests, Init)
    {
        SBSearchDatabasePostgreSQL sbsdb("dbname=sbsearch_test");
        sbsdb.close();
    }
}
