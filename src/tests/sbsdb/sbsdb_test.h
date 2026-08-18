#ifndef SBSDB_SBSDB_TEST_H_
#define SBSDB_SBSDB_TEST_H_

#include <gtest/gtest.h>

#include "sbsdb/postgresql.h"

using namespace sbsearch;
using namespace sbsearch::sbsdb;

namespace sbsearch::sbsdb::testing
{
    class SBSearchDatabaseTest : public ::testing::Test
    {
    protected:
        void SetUp() override
        {
            db.execute(R"(
                DROP TABLE IF EXISTS observation_terms_index;
                DROP TABLE IF EXISTS configuration;
                DROP TABLE IF EXISTS found;
                DROP TABLE IF EXISTS ephemerides;
                DROP TABLE IF EXISTS observatories;
                DROP TABLE IF EXISTS moving_targets;
                DROP TABLE IF EXISTS observations;
            )");
            db.setup_tables();
        }

        Postgresql db{"dbname=sbsearch_test"};
    };
}

#endif // SBSDB_SBSDB_TEST_H_