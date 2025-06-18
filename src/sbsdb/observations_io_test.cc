#include <string>
#include <vector>
#include <gtest/gtest.h>

#include "../exceptions.h"
#include "../observation.h"
#include "add.h"
#include "count.h"
#include "get.h"
#include "remove.h"
#include "sbsdb_test.h"
#include "update.h"

using namespace sbsearch;
using namespace sbsearch::sbsdb;
using std::string;
using std::vector;

namespace sbsearch::sbsdb::testing
{
    // SBSearchDatabaseTest fixture defined in sbsdb_test.cc

    TEST_F(SBSearchDatabaseTest, ObservationsIO)
    {
        // nothing in the database yet
        EXPECT_EQ(count::observations(&db, 0, 10), 0);

        Observations obs{{"test source", "X05", "product", 0, 1, "0:0, 0:1, 1:1"}};

        // verify that observation_id is not yet defined
        EXPECT_NO_THROW(verify::observations(obs, false, false));

        // but required terms are not yet defined
        EXPECT_THROW(verify::observations(obs, false, true), ObservationError);

        // update terms, but center still missing
        obs[0].terms(vector<string>({"asdf", "fdsa"}));
        EXPECT_THROW(verify::observations(obs, false, true), ObservationError);

        // update center, should be good now
        obs[0].center("1");
        EXPECT_NO_THROW(verify::observations(obs, false, true));

        // add observation, and check updated observation_id
        add::observations(&db, obs);
        EXPECT_TRUE(obs[0].observation_id());
        EXPECT_EQ(count::observations(&db, 0, 10), 1);

        // get that observation from the database and check that it matches
        Observations retrieved = get::observations(&db, obs.observation_ids());
        EXPECT_TRUE(retrieved[0] == obs[0]);

        // before updating, verify that observation_ids and terms are defined
        EXPECT_NO_THROW(verify::observations(obs, true, true));

        // edit the observation and update
        obs[0].terms(vector<string>({"a", "b", "c"}));
        update::observations(&db, {obs});
        retrieved = get::observations(&db, {obs[0].observation_id()});
        EXPECT_EQ(retrieved[0].terms(), vector<string>({"a", "b", "c"}));
        EXPECT_EQ(count::observations(&db, 0, 10), 1);

        // add a different source
        obs = Observations{{"another test source", "T05", "d", 4, 5, "0:0, 0:1, 1:1", {"d", "e", "f"}, {}, "e"}};
        add::observations(&db, obs);
        EXPECT_EQ(count::observations(&db, 0, 3), 1);
        EXPECT_EQ(count::observations(&db, 3, 10), 1);
        EXPECT_EQ(count::observations(&db, 0, 10), 2);
        EXPECT_EQ(count::observations(&db, "test source", 0, 10), 1);
        EXPECT_EQ(count::observations(&db, "another test source", 0, 10), 1);

        // remove an observation and try to get it from the database.
        auto ids = obs.observation_ids();
        remove::observations(&db, obs);
        EXPECT_FALSE(obs[0].observation_id());
        EXPECT_THROW(get::observations(&db, ids), std::runtime_error);
    }

    TEST_F(SBSearchDatabaseTest, GetAllObservations)
    {
        Observations observations({
            Observation("test source 1", "X05", "product1", 0, 1, "0:0, 0:1, 1:1", {"a", "b", "c"}, {}, "b"),
            Observation("test source 2", "568", "product2", 1, 2, "0:1, 0:2, 1:2", {"b", "c", "d"}, {}, "c"),
            Observation("test source 1", "X05", "product3", 2, 3, "0:2, 0:3, 1:3", {"c", "d", "e"}, {}, "d"),
            Observation("test source 2", "568", "product4", 3, 4, "0:3, 0:4, 1:4", {"d", "e", "f"}, {}, "e"),
        });
        add::observations(&db, observations);

        auto retrieved = db.get_all_observations("observations");
        EXPECT_EQ(retrieved.size(), 4);
        EXPECT_EQ(retrieved[0].terms(), vector<string>({"a", "b", "c"}));
    }

    TEST_F(SBSearchDatabaseTest, AllObservationsFOV)
    {
        Observations observations({
            Observation("test source 1", "X05", "product1", 0, 1, "0:0, 0:1, 1:1", {"a", "b", "c"}, {}, "b"),
            Observation("test source 2", "568", "product2", 1, 2, "0:1, 0:2, 1:2", {"b", "c", "d"}, {}, "c"),
            Observation("test source 1", "X05", "product3", 2, 3, "0:2, 0:3, 1:3", {"c", "d", "e"}, {}, "d"),
            Observation("test source 2", "568", "product4", 3, 4, "0:3, 0:4, 1:4", {"d", "e", "f"}, {}, "e"),
        });
        add::observations(&db, observations);

        auto rows = get::all_observations_fov(&db, 2, 0);
        EXPECT_EQ(rows.size(), 2);
        EXPECT_EQ(rows[0].first, observations[0].observation_id());
        EXPECT_EQ(rows[0].second, "0:0, 0:1, 1:1");
        EXPECT_EQ(rows[1].first, observations[1].observation_id());
        EXPECT_EQ(rows[1].second, "0:1, 0:2, 1:2");

        rows = get::all_observations_fov(&db, 2, 2);
        EXPECT_EQ(rows.size(), 2);
        EXPECT_EQ(rows[0].first, observations[2].observation_id());
        EXPECT_EQ(rows[0].second, "0:2, 0:3, 1:3");
        EXPECT_EQ(rows[1].first, observations[3].observation_id());
        EXPECT_EQ(rows[1].second, "0:3, 0:4, 1:4");

        rows = get::all_observations_fov(&db, 4, 1);
        EXPECT_EQ(rows.size(), 3);
        EXPECT_EQ(rows[0].first, observations[1].observation_id());
        EXPECT_EQ(rows[0].second, "0:1, 0:2, 1:2");
        EXPECT_EQ(rows[1].first, observations[2].observation_id());
        EXPECT_EQ(rows[1].second, "0:2, 0:3, 1:3");
        EXPECT_EQ(rows[2].first, observations[3].observation_id());
        EXPECT_EQ(rows[2].second, "0:3, 0:4, 1:4");

        rows = get::all_observations_fov(&db, 4, 4);
        EXPECT_EQ(rows.size(), 0);
    }

    TEST_F(SBSearchDatabaseTest, ObservationDateRange)
    {
        Observations observations({
            Observation("test source 1", "X05", "product1", 0, 1, "0:0, 0:1, 1:1"),
            Observation("test source 2", "568", "product2", 1, 2, "0:0, 0:1, 1:1"),
            Observation("test source 1", "X05", "product3", 2, 3, "0:0, 0:1, 1:1"),
            Observation("test source 2", "568", "product4", 3, 4, "0:0, 0:1, 1:1"),
        });
        for (int i = 0; i < 4; i++)
        {
            observations[i].terms("asdf fdsa");
            observations[i].center("1");
        }

        add::observations(&db, observations);

        auto drange = get::observations_date_range(&db);
        EXPECT_EQ(drange.first, 0);
        EXPECT_EQ(drange.second, 4);

        drange = get::observations_date_range(&db, "test source 1");
        EXPECT_EQ(drange.first, 0);
        EXPECT_EQ(drange.second, 3);

        drange = get::observations_date_range(&db, "test source 2");
        EXPECT_EQ(drange.first.value(), 1);
        EXPECT_EQ(drange.second.value(), 4);

        // null for no observations
        drange = get::observations_date_range(&db, "test source 3");
        EXPECT_FALSE(drange.first);
        EXPECT_FALSE(drange.second);
    }

}