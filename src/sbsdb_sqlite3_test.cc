#include "config.h"

#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <gtest/gtest.h>

#include "ephemeris.h"
#include "indexer.h"
#include "moving_target.h"
#include "observation.h"
#include "sbsdb.h"
#include "sbsdb_sqlite3.h"

// using sbsearch::Observation;
// using sbsearch::SBSearchDatabase;
// using sbsearch::SBSearchDatabaseSqlite3;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace sbsearch
{
    namespace testing
    {
        TEST(SBSearchDatabaseSqlite3Tests, SBSearchDatabaseSqlite3Init)
        {
            // open a temporary private file
            SBSearchDatabaseSqlite3 sbsdb("");

            // error if closed
            sbsdb.close();
            EXPECT_THROW(sbsdb.setup_tables(), std::runtime_error);

            // try to open the root directory as a database file
            EXPECT_THROW(SBSearchDatabaseSqlite3("/"), std::runtime_error);
        }

        TEST(SBSearchDatabaseSqlite3Tests, SBSearchDatabaseSqlite3SetupTables)
        {
            SBSearchDatabaseSqlite3 sbsdb(":memory:");
            Observation obs("test source", "product", 0, 1, "0:0, 0:1, 1:1");
            Indexer indexer;
            obs.terms(indexer.index_terms(obs));

            // tables are not yet setup
            EXPECT_THROW(sbsdb.add_observation(obs), std::runtime_error);

            // set them up
            sbsdb.setup_tables();
            EXPECT_NO_THROW(sbsdb.add_observation(obs));
        }

        TEST(SBSearchDatabaseSqlite3Tests, SBSearchDatabaseSqlite3DropCreateObservationsIndices)
        {
            SBSearchDatabaseSqlite3 sbsdb(":memory:");
            sbsdb.setup_tables();
            EXPECT_EQ(*sbsdb.get_int("SELECT COUNT(*) FROM sqlite_master WHERE type='index' and name='idx_observations_mjd_start';"), 1);
            EXPECT_EQ(*sbsdb.get_int("SELECT COUNT(*) FROM sqlite_master WHERE type='index' and name='idx_observations_mjd_stop';"), 1);
            sbsdb.drop_observations_indices();
            EXPECT_EQ(*sbsdb.get_int("SELECT COUNT(*) FROM sqlite_master WHERE type='index' and name='idx_observations_mjd_start';"), 0);
            EXPECT_EQ(*sbsdb.get_int("SELECT COUNT(*) FROM sqlite_master WHERE type='index' and name='idx_observations_mjd_stop';"), 0);
            sbsdb.create_observations_indices();
            EXPECT_EQ(*sbsdb.get_int("SELECT COUNT(*) FROM sqlite_master WHERE type='index' and name='idx_observations_mjd_start';"), 1);
            EXPECT_EQ(*sbsdb.get_int("SELECT COUNT(*) FROM sqlite_master WHERE type='index' and name='idx_observations_mjd_stop';"), 1);
        }

        TEST(SBSearchDatabaseSqlite3Tests, SBSearchDatabaseSqlite3GetInt64)
        {
            SBSearchDatabaseSqlite3 sbsdb(":memory:");
            Indexer indexer;
            Observation obs("test source", "product", 1, 2, "0:0, 0:1, 1:1");
            obs.terms(indexer.index_terms(obs));

            sbsdb.setup_tables();
            sbsdb.add_observation(obs);
            int64 *value = sbsdb.get_int64("SELECT observation_id FROM observations LIMIT 1");
            EXPECT_EQ(*value, 1);

            value = sbsdb.get_int64("SELECT observation_id FROM observations WHERE source = 'invalid source'");
            EXPECT_TRUE(value == nullptr);

            // try to get a value from a table that does not exist
            EXPECT_THROW(sbsdb.get_int64("SELECT observation_id FROM invalid_table LIMIT 1"),
                         std::runtime_error);
        }

        TEST(SBSearchDatabaseSqlite3Tests, SBSearchDatabaseSqlite3GetInt)
        {
            SBSearchDatabaseSqlite3 sbsdb(":memory:");
            Indexer indexer;
            Observation obs("test source", "product", 1, 2, "0:0, 0:1, 1:1");
            obs.terms(indexer.index_terms(obs));

            sbsdb.setup_tables();
            sbsdb.add_observation(obs);
            int *value = sbsdb.get_int("SELECT observation_id FROM observations LIMIT 1");
            EXPECT_EQ(*value, 1);

            value = sbsdb.get_int("SELECT observation_id FROM observations WHERE source = 'invalid source'");
            EXPECT_TRUE(value == nullptr);

            // try to get a value from a table that does not exist
            EXPECT_THROW(sbsdb.get_int("SELECT observation_id FROM invalid_table LIMIT 1"),
                         std::runtime_error);
        }

        TEST(SBSearchDatabaseSqlite3Tests, SBSearchDatabaseSqlite3GetDouble)
        {
            SBSearchDatabaseSqlite3 sbsdb(":memory:");
            Indexer indexer;
            Observation obs("test source", "product", 1, 2, "0:0, 0:1, 1:1");
            obs.terms(indexer.index_terms(obs));

            sbsdb.setup_tables();
            sbsdb.add_observation(obs);
            double *value = sbsdb.get_double("SELECT mjd_start FROM observations LIMIT 1");
            EXPECT_EQ(*value, 1);

            value = sbsdb.get_double("SELECT mjd_start FROM observations WHERE source = 'invalid source'");
            EXPECT_TRUE(value == nullptr);

            // try to get a value from a table that does not exist
            EXPECT_THROW(sbsdb.get_double("SELECT mjd_start FROM invalid_table LIMIT 1"),
                         std::runtime_error);
        }

        TEST(SBSearchDatabaseSqlite3Tests, SBSearchDatabaseSqlite3GetString)
        {
            SBSearchDatabaseSqlite3 sbsdb(":memory:");
            Indexer indexer;
            Observation obs("test source", "product", 1, 2, "0:0, 0:1, 1:1");
            obs.terms(indexer.index_terms(obs));

            sbsdb.setup_tables();

            std::string *s = sbsdb.get_string("SELECT value FROM configuration WHERE parameter='max_spatial_cells'");
            EXPECT_EQ(*s, "8");

            s = sbsdb.get_string("SELECT value FROM configuration WHERE parameter='invalid parameter'");
            EXPECT_TRUE(s == nullptr);

            // try to get a value from a table that does not exist
            EXPECT_THROW(sbsdb.get_string("SELECT whatever FROM invalid_table LIMIT 1"),
                         std::runtime_error);
        }

        TEST(SBSearchDatabaseSqlite3Tests, SBSearchDatabaseSqlite3DateRange)
        {
            SBSearchDatabaseSqlite3 sbsdb(":memory:");
            sbsdb.setup_tables();

            std::vector<Observation> observations = {
                Observation("test source 1", "product1", 0, 1, "0:0, 0:1, 1:1"),
                Observation("test source 2", "product2", 1, 2, "0:0, 0:1, 1:1"),
                Observation("test source 1", "product3", 2, 3, "0:0, 0:1, 1:1"),
                Observation("test source 2", "product4", 3, 4, "0:0, 0:1, 1:1"),
            };
            for (int i = 0; i < 4; i++)
                observations[i].terms("asdf fdsa");

            sbsdb.add_observations(observations);

            auto drange = sbsdb.date_range();
            EXPECT_EQ(*drange.first, 0);
            EXPECT_EQ(*drange.second, 4);

            drange = sbsdb.date_range("test source 1");
            EXPECT_EQ(*drange.first, 0);
            EXPECT_EQ(*drange.second, 3);

            drange = sbsdb.date_range("test source 2");
            EXPECT_EQ(*drange.first, 1);
            EXPECT_EQ(*drange.second, 4);

            // null pointer for no observations
            drange = sbsdb.date_range("test source 3");
            EXPECT_EQ(drange.first, nullptr);
            EXPECT_EQ(drange.second, nullptr);
        }

        TEST(SBSearchDatabaseSqlite3Tests, SBSearchDatabaseSqlite3AddGetMovingTarget)
        {
            SBSearchDatabaseSqlite3 sbsdb(":memory:");
            sbsdb.setup_tables();
            MovingTarget encke("2P");
            MovingTarget ceres("1");

            // add to the database, expect an updated object_id
            EXPECT_EQ(encke.object_id(), UNDEF_OBJECT_ID);
            sbsdb.add_moving_target(encke);
            EXPECT_EQ(encke.object_id(), 1);

            // add another, it should be 2
            sbsdb.add_moving_target(ceres);
            EXPECT_EQ(ceres.object_id(), 2);

            // get them from the database
            MovingTarget test;
            test = sbsdb.get_moving_target(1);
            EXPECT_EQ(test.designation(), encke.designation());
            EXPECT_EQ(test.object_id(), encke.object_id());
            EXPECT_EQ(test.alternate_names(), encke.alternate_names());

            test = sbsdb.get_moving_target("2P");
            EXPECT_EQ(test.designation(), encke.designation());
            EXPECT_EQ(test.object_id(), encke.object_id());
            EXPECT_EQ(test.alternate_names(), encke.alternate_names());

            test = sbsdb.get_moving_target("1");
            EXPECT_EQ(test.designation(), ceres.designation());
            EXPECT_EQ(test.object_id(), ceres.object_id());
            EXPECT_EQ(test.alternate_names(), ceres.alternate_names());

            test = sbsdb.get_moving_target(2);
            EXPECT_EQ(test.designation(), ceres.designation());
            EXPECT_EQ(test.object_id(), ceres.object_id());
            EXPECT_EQ(test.alternate_names(), ceres.alternate_names());

            // add an alternate name and update encke
            encke.add_name("2P/Encke");
            sbsdb.update_moving_target(encke);
            test = sbsdb.get_moving_target("2P");
            EXPECT_EQ(test.designation(), encke.designation());
            EXPECT_EQ(test.object_id(), encke.object_id());
            EXPECT_EQ(test.alternate_names(), encke.alternate_names());

            // try getting encke via alt name
            test = sbsdb.get_moving_target("2P/Encke");
            EXPECT_EQ(test.designation(), encke.designation());
            EXPECT_EQ(test.object_id(), encke.object_id());
            EXPECT_EQ(test.alternate_names(), encke.alternate_names());

            // add a few names to ceres
            vector<string> names{"(1) Ceres", "Ceres", "A801 AA"};
            ceres.add_names(names.begin(), names.end());
            sbsdb.update_moving_target(ceres);
            test = sbsdb.get_moving_target("A801 AA");
            EXPECT_EQ(test.designation(), ceres.designation());
            EXPECT_EQ(test.object_id(), ceres.object_id());
            EXPECT_EQ(test.alternate_names(), ceres.alternate_names());
            EXPECT_EQ(test.alternate_names().size(), 3);

            // Try to add an object that already exists.
            MovingTarget halley("1P");
            halley.object_id(1);
            EXPECT_THROW(sbsdb.add_moving_target(halley), MovingTargetError);

            // Try to remove an object that does not exist.
            MovingTarget new_comet("1000P");
            new_comet.object_id(9123);
            EXPECT_THROW(sbsdb.remove_moving_target(new_comet), MovingTargetError);

            // Get an object that does not exist
            EXPECT_THROW(sbsdb.get_moving_target(1000), MovingTargetError);
            EXPECT_THROW(sbsdb.get_moving_target("asdf"), MovingTargetError);
        }

        TEST(SBSearchDatabaseSqlite3Tests, SBSearchDatabaseSqlite3AddGetObservation)
        {
            SBSearchDatabaseSqlite3 sbsdb(":memory:");
            sbsdb.setup_tables();

            Observation obs("test source", "product", 0, 1, "0:0, 0:1, 1:1");
            // observation_id is not yet defined
            EXPECT_EQ(obs.observation_id(), UNDEFINED_OBSID);

            // terms are not yet defined
            EXPECT_THROW(sbsdb.add_observation(obs), std::runtime_error);

            // update terms, add observation, now observation_id should be updated
            obs.terms(vector<string>{"asdf", "fdsa"});
            sbsdb.add_observation(obs);
            EXPECT_NE(obs.observation_id(), UNDEFINED_OBSID);

            Observation retrieved = sbsdb.get_observation(obs.observation_id());
            EXPECT_TRUE(retrieved.is_equal(obs));

            // edit the observation and update
            obs.terms(vector<string>{"a", "b", "c"});
            sbsdb.add_observation(obs);
            retrieved = sbsdb.get_observation(obs.observation_id());
            EXPECT_EQ(retrieved.terms(), "a b c");

            // try to get an observation that does not exist
            EXPECT_THROW(sbsdb.get_observation(-1), std::runtime_error);
        }

        TEST(SBSearchDatabaseSqlite3Tests, SBSearchDatabaseSqlite3FindObservations)
        {

            SBSearchDatabaseSqlite3 sbsdb(":memory:");
            sbsdb.setup_tables();

            Observation obs("test source", "a", 0, 1, "0:0, 0:1, 1:1", "a b c");
            sbsdb.add_observation(obs);

            obs = Observation("test source", "b", 0, 1, "0:0, 0:1, 1:1", "b c d");
            sbsdb.add_observation(obs);

            obs = Observation("test source", "c", 0, 1, "0:0, 0:1, 1:1", "c d e");
            sbsdb.add_observation(obs);

            obs = Observation("test source", "d", 0, 1, "0:0, 0:1, 1:1", "d e f");
            sbsdb.add_observation(obs);

            // find observations matching term a
            vector<Observation> matches;
            matches = sbsdb.find_observations(vector<string>{"a"});
            EXPECT_EQ(matches.size(), 1);

            // a or f
            matches = sbsdb.find_observations(vector<string>{"a", "f"});
            EXPECT_EQ(matches.size(), 2);

            // c or f
            matches = sbsdb.find_observations(vector<string>{"c", "f"});
            EXPECT_EQ(matches.size(), 4);

            // g
            matches = sbsdb.find_observations(vector<string>{"g"});
            EXPECT_EQ(matches.size(), 0);
        }

        TEST(SBSearchDatabaseSqlite3Tests, SBSearchDatabaseSqlite3CheckSQL)
        {
            SBSearchDatabaseSqlite3 sbsdb(":memory:");
        }

        TEST(SBSearchDatabaseSqlite3Tests, SBSearchDatabaseSqlite3ErrorIfClosed)
        {
            SBSearchDatabaseSqlite3 sbsdb(":memory:");
            sbsdb.close();
        }
    }
}
