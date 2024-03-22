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

using std::cerr;
using std::cout;
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
            Observation obs("test source", "X05", "product", 0, 1, "0:0, 0:1, 1:1");
            Indexer indexer;
            obs.terms(indexer.index_terms(obs));
            EXPECT_NO_THROW(sbsdb.add_observation(obs));
            EXPECT_NO_THROW(sbsdb.setup_tables());
        }

        TEST(SBSearchDatabaseSqlite3Tests, SBSearchDatabaseSqlite3DropCreateObservationsIndices)
        {
            SBSearchDatabaseSqlite3 sbsdb(":memory:");
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
            Observation obs("test source", "X05", "product", 1, 2, "0:0, 0:1, 1:1");
            obs.terms(indexer.index_terms(obs));

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
            Observation obs("test source", "X05", "product", 1, 2, "0:0, 0:1, 1:1");
            obs.terms(indexer.index_terms(obs));

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
            Observation obs("test source", "X05", "product", 1, 2, "0:0, 0:1, 1:1");
            obs.terms(indexer.index_terms(obs));

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
            Observation obs("test source", "X05", "product", 1, 2, "0:0, 0:1, 1:1");
            obs.terms(indexer.index_terms(obs));

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

            Observations observations = {
                Observation("test source 1", "X05", "product1", 0, 1, "0:0, 0:1, 1:1"),
                Observation("test source 2", "568", "product2", 1, 2, "0:0, 0:1, 1:1"),
                Observation("test source 1", "X05", "product3", 2, 3, "0:0, 0:1, 1:1"),
                Observation("test source 2", "568", "product4", 3, 4, "0:0, 0:1, 1:1"),
            };
            for (int i = 0; i < 4; i++)
                observations[i].terms("asdf fdsa");

            sbsdb.add_observations(observations);

            auto drange = sbsdb.observation_date_range();
            EXPECT_EQ(*drange.first, 0);
            EXPECT_EQ(*drange.second, 4);

            drange = sbsdb.observation_date_range("test source 1");
            EXPECT_EQ(*drange.first, 0);
            EXPECT_EQ(*drange.second, 3);

            drange = sbsdb.observation_date_range("test source 2");
            EXPECT_EQ(*drange.first, 1);
            EXPECT_EQ(*drange.second, 4);

            // null pointer for no observations
            drange = sbsdb.observation_date_range("test source 3");
            EXPECT_EQ(drange.first, nullptr);
            EXPECT_EQ(drange.second, nullptr);
        }

        TEST(SBSearchDatabaseSqlite3Tests, SBSearchDatabaseSqlite3AddGetMovingTarget)
        {
            SBSearchDatabaseSqlite3 sbsdb(":memory:");
            MovingTarget encke("2P");
            MovingTarget ceres("1");
            MovingTarget mercury("1", false);

            // add to the database, expect an updated moving_target_id
            EXPECT_EQ(encke.moving_target_id(), UNDEF_MOVING_TARGET_ID);
            sbsdb.add_moving_target(encke);
            EXPECT_EQ(encke.moving_target_id(), 1);

            // add another, it should be 2
            sbsdb.add_moving_target(ceres);
            EXPECT_EQ(ceres.moving_target_id(), 2);

            // and this should be 3
            sbsdb.add_moving_target(mercury);
            EXPECT_EQ(mercury.moving_target_id(), 3);

            // get them from the database
            MovingTarget test;
            test = sbsdb.get_moving_target(1);
            EXPECT_EQ(test.designation(), encke.designation());
            EXPECT_EQ(test.moving_target_id(), encke.moving_target_id());
            EXPECT_EQ(test.alternate_names(), encke.alternate_names());
            EXPECT_EQ(test.small_body(), encke.small_body());

            test = sbsdb.get_moving_target("2P");
            EXPECT_EQ(test.designation(), encke.designation());
            EXPECT_EQ(test.moving_target_id(), encke.moving_target_id());
            EXPECT_EQ(test.alternate_names(), encke.alternate_names());
            EXPECT_EQ(test.small_body(), encke.small_body());

            test = sbsdb.get_moving_target("1");
            EXPECT_EQ(test.designation(), ceres.designation());
            EXPECT_EQ(test.moving_target_id(), ceres.moving_target_id());
            EXPECT_EQ(test.alternate_names(), ceres.alternate_names());
            EXPECT_EQ(test.small_body(), ceres.small_body());

            test = sbsdb.get_moving_target(2);
            EXPECT_EQ(test.designation(), ceres.designation());
            EXPECT_EQ(test.moving_target_id(), ceres.moving_target_id());
            EXPECT_EQ(test.alternate_names(), ceres.alternate_names());
            EXPECT_EQ(test.small_body(), ceres.small_body());

            test = sbsdb.get_moving_target("1", false);
            EXPECT_EQ(test.designation(), mercury.designation());
            EXPECT_EQ(test.moving_target_id(), mercury.moving_target_id());
            EXPECT_EQ(test.alternate_names(), mercury.alternate_names());
            EXPECT_EQ(test.small_body(), mercury.small_body());

            test = sbsdb.get_moving_target(3);
            EXPECT_EQ(test.designation(), mercury.designation());
            EXPECT_EQ(test.moving_target_id(), mercury.moving_target_id());
            EXPECT_EQ(test.alternate_names(), mercury.alternate_names());
            EXPECT_EQ(test.small_body(), mercury.small_body());

            // add an alternate name and update encke
            encke.add_name("2P/Encke");
            sbsdb.update_moving_target(encke);
            test = sbsdb.get_moving_target("2P");
            EXPECT_EQ(test.designation(), encke.designation());
            EXPECT_EQ(test.moving_target_id(), encke.moving_target_id());
            EXPECT_EQ(test.alternate_names(), encke.alternate_names());

            // try getting Encke via alt name
            test = sbsdb.get_moving_target("2P/Encke");
            EXPECT_EQ(test.designation(), encke.designation());
            EXPECT_EQ(test.moving_target_id(), encke.moving_target_id());
            EXPECT_EQ(test.alternate_names(), encke.alternate_names());

            // add a few names to Ceres
            vector<string> names{"(1) Ceres", "Ceres", "A801 AA"};
            ceres.add_names(names.begin(), names.end());
            sbsdb.update_moving_target(ceres);
            test = sbsdb.get_moving_target("A801 AA");
            EXPECT_EQ(test.designation(), ceres.designation());
            EXPECT_EQ(test.moving_target_id(), ceres.moving_target_id());
            EXPECT_EQ(test.alternate_names(), ceres.alternate_names());
            EXPECT_EQ(test.alternate_names().size(), 3);

            // and Mercury
            names = vector<string>{"Hermes", "Nabu"};
            mercury.add_names(names.begin(), names.end());
            sbsdb.update_moving_target(mercury);
            test = sbsdb.get_moving_target("Nabu", false);
            EXPECT_EQ(test.designation(), mercury.designation());
            EXPECT_EQ(test.moving_target_id(), mercury.moving_target_id());
            EXPECT_EQ(test.alternate_names(), mercury.alternate_names());
            EXPECT_EQ(test.alternate_names().size(), 2);

            // Try to add an object that already exists.
            MovingTarget halley("1P");
            halley.moving_target_id(1);
            EXPECT_THROW(sbsdb.add_moving_target(halley), MovingTargetError);

            MovingTarget duplicate_ceres("1");
            EXPECT_THROW(sbsdb.add_moving_target(duplicate_ceres), MovingTargetError);

            // Try to remove an object that does not exist.
            MovingTarget new_comet("1000P");
            new_comet.moving_target_id(9123);
            EXPECT_THROW(sbsdb.remove_moving_target(new_comet), MovingTargetError);

            // Get objects that do not exist
            EXPECT_THROW(sbsdb.get_moving_target(1000), MovingTargetError);
            EXPECT_EQ(sbsdb.get_moving_target("asdf"), MovingTarget("asdf"));
            EXPECT_EQ(sbsdb.get_moving_target("Nabu", true).moving_target_id(), UNDEF_MOVING_TARGET_ID);
            EXPECT_EQ(sbsdb.get_moving_target("Ceres", false).moving_target_id(), UNDEF_MOVING_TARGET_ID);
        }

        TEST(SBSearchDatabaseSqlite3Tests, SBSearchDatabaseSqlite3AddGetObservatory)
        {
            SBSearchDatabaseSqlite3 sbsdb(":memory:");

            const Observatory ztf{243.14022, 0.836325, +0.546877};
            const Observatory ldt{248.57749, 0.822887, 0.566916};
            const Observatory maunakea{204.5278, 0.94171, +0.33725};
            const Observatory paranal{289.59569, 0.909943, -0.414336};

            sbsdb.add_observatory("I41", ztf);
            sbsdb.add_observatory("G37", ldt);
            sbsdb.add_observatory("568", maunakea);
            sbsdb.add_observatory("309", paranal);

            Observatory obs = sbsdb.get_observatory("I41");
            EXPECT_EQ(obs, ztf);

            obs = sbsdb.get_observatory("G37");
            EXPECT_EQ(obs, ldt);

            obs = sbsdb.get_observatory("568");
            EXPECT_EQ(obs, maunakea);

            obs = sbsdb.get_observatory("309");
            EXPECT_EQ(obs, paranal);

            Observatories observatories = sbsdb.get_observatories();
            EXPECT_EQ(observatories["I41"], ztf);
            EXPECT_EQ(observatories["G37"], ldt);
            EXPECT_EQ(observatories["568"], maunakea);
            EXPECT_EQ(observatories["309"], paranal);

            sbsdb.remove_observatory("G37");
            EXPECT_THROW(sbsdb.get_observatory("G37"), ObservatoryError);

            EXPECT_THROW(sbsdb.add_observatory("I41", ztf), ObservatoryError);
        }

        TEST(SBSearchDatabaseSqlite3Tests, SBSearchDatabaseSqlite3AddGetEphemeris)
        {
            SBSearchDatabaseSqlite3 sbsdb(":memory:");

            MovingTarget encke{"2P"};
            Ephemeris eph{encke,
                          {{0, 10, 1, 0, 1, 0.1, 90, 0, 1, 180, 0, 0, 0, 10, -1},
                           {1, 11, 2, 0, 5, 0.5, 90, 1, 0, 0, 180, 30, 0, 20, 5},
                           {2, 12, 3, 0, 10, 1.0, 90, 2, 1, 90, 80, 90, 0, 30, 10}}};

            // The target is not in the database, so we expect an error
            EXPECT_THROW(sbsdb.add_ephemeris(eph), MovingTargetError);

            // Add the target, verify that the id was updated
            sbsdb.add_moving_target(encke);
            EXPECT_NE(encke.moving_target_id(), eph.target().moving_target_id());

            // Fix the target, and then we can add the ephemeris data
            eph.target(encke);
            sbsdb.add_ephemeris(eph);

            // Get the data back
            Ephemeris test;
            test = sbsdb.get_ephemeris(eph.target());
            EXPECT_EQ(test, eph);

            // Get a subset of data
            test = sbsdb.get_ephemeris(eph.target(), 0.5, 1.5);
            EXPECT_EQ(test, eph[1]);

            // This target does not match database copy
            MovingTarget wrong_id{"1P", eph.target().moving_target_id()};
            Ephemeris other{wrong_id, eph.data()};
            EXPECT_THROW(sbsdb.add_ephemeris(other), MovingTargetError);

            // Remove some data
            sbsdb.remove_ephemeris(eph.target(), 1.5, 10);
            test = sbsdb.get_ephemeris(eph.target());
            EXPECT_NE(test, eph);
            EXPECT_EQ(test, eph.slice(0, 2));

            // Remove all
            sbsdb.remove_ephemeris(eph.target());
            test = sbsdb.get_ephemeris(eph.target());
            EXPECT_EQ(test.num_vertices(), 0);
        }

        TEST(SBSearchDatabaseSqlite3Tests, SBSearchDatabaseSqlite3AddGetObservation)
        {
            SBSearchDatabaseSqlite3 sbsdb(":memory:");

            Observation obs("test source", "X05", "product", 0, 1, "0:0, 0:1, 1:1");
            // observation_id is not yet defined
            EXPECT_EQ(obs.observation_id(), UNDEFINED_OBSID);

            // terms are not yet defined
            EXPECT_THROW(sbsdb.add_observation(obs), std::runtime_error);

            // update terms, add observation, now observation_id should be updated
            obs.terms(vector<string>{"asdf", "fdsa"});
            sbsdb.add_observation(obs);
            EXPECT_NE(obs.observation_id(), UNDEFINED_OBSID);

            Observation retrieved = sbsdb.get_observation(obs.observation_id());
            EXPECT_TRUE(retrieved == obs);

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

            Observation obs("test source", "X05", "a", 0, 1, "0:0, 0:1, 1:1", "a b c");
            sbsdb.add_observation(obs);

            obs = Observation("test source", "X05", "b", 1, 2, "0:0, 0:1, 1:1", "b c d");
            sbsdb.add_observation(obs);

            obs = Observation("test source", "X05", "c", 2, 3, "0:0, 0:1, 1:1", "c d e");
            sbsdb.add_observation(obs);

            obs = Observation("another test source", "T05", "d", 4, 5, "0:0, 0:1, 1:1", "d e f");
            sbsdb.add_observation(obs);

            // find observations matching term a
            Observations matches;
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

            // test observation time limits
            // start
            matches = sbsdb.find_observations({"e"}, {.mjd_start = 2});
            EXPECT_EQ(matches.size(), 2);

            matches = sbsdb.find_observations({"e"}, {.mjd_start = 3.5});
            EXPECT_EQ(matches.size(), 1);

            // stop
            matches = sbsdb.find_observations({"e"}, {.mjd_stop = 1});
            EXPECT_EQ(matches.size(), 0);

            matches = sbsdb.find_observations({"e"}, {.mjd_stop = 3});
            EXPECT_EQ(matches.size(), 1);

            matches = sbsdb.find_observations({"e"}, {.mjd_stop = 5});
            EXPECT_EQ(matches.size(), 2);

            // start-stop
            matches = sbsdb.find_observations({"e"}, {.mjd_start = 2, .mjd_stop = 2.5});
            EXPECT_EQ(matches.size(), 0);

            matches = sbsdb.find_observations({"e"}, {.mjd_start = 2, .mjd_stop = 3});
            EXPECT_EQ(matches.size(), 1);

            matches = sbsdb.find_observations({"e"}, {.mjd_start = 2.5, .mjd_stop = 4.5});
            EXPECT_EQ(matches.size(), 0);

            matches = sbsdb.find_observations({"e"}, {.mjd_start = 3, .mjd_stop = 5});
            EXPECT_EQ(matches.size(), 1);

            // search by source
            matches = sbsdb.find_observations({"b", "e"}, {.source = "test source"});
            EXPECT_EQ(matches.size(), 3);

            matches = sbsdb.find_observations({"b", "e"}, {.source = "another test source"});
            EXPECT_EQ(matches.size(), 1);
        }

        TEST(SBSearchDatabaseSqlite3Tests, SBSearchDatabaseSqlite3AddGetFound)
        {
            SBSearchDatabaseSqlite3 sbsdb(":memory:");

            Observations obs{
                {"test source", "X05", "a", 0, 1, "0:0, 0:1, 1:1", "a b c"},
                {"test source", "X05", "b", 1, 2, "0:0, 0:1, 1:1", "b c d"}};
            sbsdb.add_observations(obs);

            MovingTarget encke("2P");
            sbsdb.add_moving_target(encke);

            Ephemeris eph{encke,
                          {{0, 10, 1, 0, 1, 0.1, 90, 0, 1, 180, 0, 0, 0, 10, -1},
                           {1, 11, 2, 0, 5, 0.5, 90, 1, 0, 0, 180, 30, 0, 20, 5},
                           {2, 12, 3, 0, 10, 1.0, 90, 2, 1, 90, 80, 90, 0, 30, 10}}};

            // these may not make sense, but it doesn't matter
            vector<Found> founds = {
                {obs[0], eph[0]},
                {obs[1], eph[1]}};
            sbsdb.add_founds(founds);

            founds = sbsdb.get_found(obs[0]);
            EXPECT_EQ(founds.size(), 1);
            EXPECT_EQ(founds[0], Found(obs[0], eph[0]));

            founds = sbsdb.get_found(obs[1]);
            EXPECT_EQ(founds.size(), 1);
            EXPECT_EQ(founds[0], Found(obs[1], eph[1]));

            founds = sbsdb.get_found(encke);
            EXPECT_EQ(founds.size(), 2);
            EXPECT_EQ(std::count(founds.begin(), founds.end(), Found(obs[0], eph[0])), 1);
            EXPECT_EQ(std::count(founds.begin(), founds.end(), Found(obs[1], eph[1])), 1);

            sbsdb.remove_founds(founds);
            founds = sbsdb.get_found(obs[0]);
            EXPECT_EQ(founds.size(), 0);
            founds = sbsdb.get_found(obs[1]);
            EXPECT_EQ(founds.size(), 0);
            founds = sbsdb.get_found(encke);
            EXPECT_EQ(founds.size(), 0);
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
