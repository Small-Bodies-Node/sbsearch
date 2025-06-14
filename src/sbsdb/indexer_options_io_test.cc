#include <gtest/gtest.h>

#include "../indexer.h"
#include "get.h"
#include "sbsdb_test.h"
#include "update.h"

using namespace sbsearch;
using namespace sbsearch::sbsdb;

namespace sbsearch::sbsdb::testing
{
    // SBSearchDatabaseTest fixture defined in sbsdb_test.cc

    TEST_F(SBSearchDatabaseTest, IndexerOptionsIO)
    {
        auto options = get::indexer_options(&db);

        // expect the defaults
        EXPECT_EQ(options.max_spatial_index_cells(), 8);
        EXPECT_EQ(options.max_spatial_level(), 12);
        EXPECT_EQ(options.min_spatial_level(), 4);

        // change everything and verify
        Indexer::Options new_options;
        new_options.max_spatial_index_cells(10);
        new_options.max_spatial_level(11);
        new_options.min_spatial_level(6);
        update::indexer_options(&db, new_options);
        EXPECT_TRUE(new_options == get::indexer_options(&db));
    }

}