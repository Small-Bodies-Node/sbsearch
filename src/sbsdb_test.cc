#include "config.h"

#include <iostream>
#include <memory>
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
#include "sbsdb_sqlite3.h"

using namespace sbsearch;
using std::cerr;
using std::cout;
using std::string;
using std::vector;

namespace testing
{
    using ::testing::TestWithParam;
    using ::testing::Values;

    typedef std::unique_ptr<SBSearchDatabase> CreateSBSearchDatabaseFunc();

    std::unique_ptr<SBSearchDatabase> CreateSqlite3Database()
    {
        return std::make_unique<SBSearchDatabaseSqlite3>(":memory:");
    }

    std::unique_ptr<SBSearchDatabase> CreatePostgreSQLDatabase()
    {
        return std::make_unique<SBSearchDatabasePostgreSQL>("dbname=sbsearch_test");
    }

    class SBSearchDatabaseTest : public TestWithParam<CreateSBSearchDatabaseFunc *>
    {
    public:
        void SetUp() override
        {
            sbsdb = (*GetParam())();
            sbsdb->execute_sql(R"(
                DROP TABLE IF EXISTS observation_terms_index;
                DROP TABLE IF EXISTS configuration;
                DROP TABLE IF EXISTS found;
                DROP TABLE IF EXISTS ephemerides;
                DROP TABLE IF EXISTS observatories;
                DROP TABLE IF EXISTS moving_targets;
                DROP TABLE IF EXISTS observations;
            )");
            sbsdb->setup_tables();
        }

    protected:
        std::unique_ptr<SBSearchDatabase> sbsdb;
    };

    TEST_P(SBSearchDatabaseTest, SetupTables)
    {
        Observation obs("test source", "X05", "product", 0, 1, "0:0, 0:1, 1:1");
        Indexer indexer;
        obs.terms(indexer.terms(Indexer::index, obs));
        EXPECT_NO_THROW(sbsdb->add_observations(obs));
        EXPECT_NO_THROW(sbsdb->setup_tables());
    }

    TEST_P(SBSearchDatabaseTest, GetInt64)
    {
        auto value = sbsdb->get_int64("SELECT 1");
        EXPECT_EQ(value, 1);

        value = sbsdb->get_int64("SELECT NULL");
        EXPECT_FALSE(value);

        // try to get a value from a table that does not exist
        EXPECT_THROW(sbsdb->get_int64("SELECT observation_id FROM invalid_table LIMIT 1"),
                     std::runtime_error);
    }

    TEST_P(SBSearchDatabaseTest, GetInt)
    {
        auto value = sbsdb->get_int("SELECT 1");
        EXPECT_EQ(value, 1);

        value = sbsdb->get_int("SELECT NULL");
        EXPECT_FALSE(value);

        // try to get a value from a table that does not exist
        EXPECT_THROW(sbsdb->get_int("SELECT observation_id FROM invalid_table LIMIT 1"),
                     std::runtime_error);
    }

    TEST_P(SBSearchDatabaseTest, GetDouble)
    {
        auto mjd_start = sbsdb->get_double("SELECT 1");
        EXPECT_EQ(mjd_start.value(), 1);

        mjd_start = sbsdb->get_double("SELECT NULL");
        EXPECT_FALSE(mjd_start);

        // try to get a value from a table that does not exist
        EXPECT_THROW(sbsdb->get_double("SELECT mjd_start FROM invalid_table LIMIT 1"),
                     std::runtime_error);
    }

    TEST_P(SBSearchDatabaseTest, GetString)
    {
        auto s = sbsdb->get_string("SELECT 'asdf'");
        EXPECT_EQ(s, "asdf");

        s = sbsdb->get_string("SELECT NULL");
        EXPECT_FALSE(s);

        // try to get a value from a table that does not exist
        EXPECT_THROW(sbsdb->get_string("SELECT whatever FROM invalid_table LIMIT 1"),
                     std::runtime_error);
    }

    TEST_P(SBSearchDatabaseTest, DateRange)
    {
        Observations observations({
            Observation("test source 1", "X05", "product1", 0, 1, "0:0, 0:1, 1:1"),
            Observation("test source 2", "568", "product2", 1, 2, "0:0, 0:1, 1:1"),
            Observation("test source 1", "X05", "product3", 2, 3, "0:0, 0:1, 1:1"),
            Observation("test source 2", "568", "product4", 3, 4, "0:0, 0:1, 1:1"),
        });
        for (int i = 0; i < 4; i++)
            observations[i].terms("asdf fdsa");

        sbsdb->add_observations(observations);

        auto drange = sbsdb->observation_date_range();
        EXPECT_EQ(drange.first, 0);
        EXPECT_EQ(drange.second, 4);

        drange = sbsdb->observation_date_range("test source 1");
        EXPECT_EQ(drange.first, 0);
        EXPECT_EQ(drange.second, 3);

        drange = sbsdb->observation_date_range("test source 2");
        EXPECT_EQ(drange.first.value(), 1);
        EXPECT_EQ(drange.second.value(), 4);

        // null pointer for no observations
        drange = sbsdb->observation_date_range("test source 3");
        EXPECT_FALSE(drange.first);
        EXPECT_FALSE(drange.second);
    }

    TEST_P(SBSearchDatabaseTest, AddGetMovingTarget)
    {
        MovingTarget encke("2P");
        MovingTarget ceres("1");
        MovingTarget mercury("1", false);

        // add to the database, expect an updated moving_target_id
        EXPECT_FALSE(encke.moving_target_id());
        sbsdb->add_moving_target(encke);
        EXPECT_EQ(encke.moving_target_id(), 1);

        // add another, it should be 2
        sbsdb->add_moving_target(ceres);
        EXPECT_EQ(ceres.moving_target_id(), 2);

        // and this should be 3
        sbsdb->add_moving_target(mercury);
        EXPECT_EQ(mercury.moving_target_id(), 3);

        // get them from the database
        MovingTarget test;
        test = sbsdb->get_moving_target(1);
        EXPECT_EQ(test.designation(), encke.designation());
        EXPECT_EQ(test.moving_target_id(), encke.moving_target_id());
        EXPECT_EQ(test.alternate_names(), encke.alternate_names());
        EXPECT_EQ(test.small_body(), encke.small_body());

        test = sbsdb->get_moving_target("2P");
        EXPECT_EQ(test.designation(), encke.designation());
        EXPECT_EQ(test.moving_target_id(), encke.moving_target_id());
        EXPECT_EQ(test.alternate_names(), encke.alternate_names());
        EXPECT_EQ(test.small_body(), encke.small_body());

        test = sbsdb->get_moving_target("1");
        EXPECT_EQ(test.designation(), ceres.designation());
        EXPECT_EQ(test.moving_target_id(), ceres.moving_target_id());
        EXPECT_EQ(test.alternate_names(), ceres.alternate_names());
        EXPECT_EQ(test.small_body(), ceres.small_body());

        test = sbsdb->get_moving_target(2);
        EXPECT_EQ(test.designation(), ceres.designation());
        EXPECT_EQ(test.moving_target_id(), ceres.moving_target_id());
        EXPECT_EQ(test.alternate_names(), ceres.alternate_names());
        EXPECT_EQ(test.small_body(), ceres.small_body());

        test = sbsdb->get_moving_target("1", false);
        EXPECT_EQ(test.designation(), mercury.designation());
        EXPECT_EQ(test.moving_target_id(), mercury.moving_target_id());
        EXPECT_EQ(test.alternate_names(), mercury.alternate_names());
        EXPECT_EQ(test.small_body(), mercury.small_body());

        test = sbsdb->get_moving_target(3);
        EXPECT_EQ(test.designation(), mercury.designation());
        EXPECT_EQ(test.moving_target_id(), mercury.moving_target_id());
        EXPECT_EQ(test.alternate_names(), mercury.alternate_names());
        EXPECT_EQ(test.small_body(), mercury.small_body());

        // add an alternate name and update encke
        encke.add_name("2P/Encke");
        sbsdb->update_moving_target(encke);
        test = sbsdb->get_moving_target("2P");
        EXPECT_EQ(test.designation(), encke.designation());
        EXPECT_EQ(test.moving_target_id(), encke.moving_target_id());
        EXPECT_EQ(test.alternate_names(), encke.alternate_names());

        // try getting Encke via alt name
        test = sbsdb->get_moving_target("2P/Encke");
        EXPECT_EQ(test.designation(), encke.designation());
        EXPECT_EQ(test.moving_target_id(), encke.moving_target_id());
        EXPECT_EQ(test.alternate_names(), encke.alternate_names());

        // add a few names to Ceres
        vector<string> names{"(1) Ceres", "Ceres", "A801 AA"};
        ceres.add_names(names.begin(), names.end());
        sbsdb->update_moving_target(ceres);
        test = sbsdb->get_moving_target("A801 AA");
        EXPECT_EQ(test.designation(), ceres.designation());
        EXPECT_EQ(test.moving_target_id(), ceres.moving_target_id());
        EXPECT_EQ(test.alternate_names(), ceres.alternate_names());
        EXPECT_EQ(test.alternate_names().size(), 3);

        // and Mercury
        names = vector<string>{"Hermes", "Nabu"};
        mercury.add_names(names.begin(), names.end());
        sbsdb->update_moving_target(mercury);
        test = sbsdb->get_moving_target("Nabu", false);
        EXPECT_EQ(test.designation(), mercury.designation());
        EXPECT_EQ(test.moving_target_id(), mercury.moving_target_id());
        EXPECT_EQ(test.alternate_names(), mercury.alternate_names());
        EXPECT_EQ(test.alternate_names().size(), 2);

        // Try to add an object that already exists.
        MovingTarget halley("1P");
        halley.moving_target_id(1);
        EXPECT_THROW(sbsdb->add_moving_target(halley), MovingTargetError);

        MovingTarget duplicate_ceres("1");
        EXPECT_THROW(sbsdb->add_moving_target(duplicate_ceres), MovingTargetError);

        // Try to remove an object that does not exist.
        MovingTarget new_comet("1000P");
        new_comet.moving_target_id(9123);
        EXPECT_THROW(sbsdb->remove_moving_target(new_comet), MovingTargetError);

        // Get objects that do not exist
        EXPECT_THROW(sbsdb->get_moving_target(1000), MovingTargetError);
        EXPECT_EQ(sbsdb->get_moving_target("asdf"), MovingTarget("asdf"));
        EXPECT_FALSE(sbsdb->get_moving_target("Nabu", true).moving_target_id());
        EXPECT_FALSE(sbsdb->get_moving_target("Ceres", false).moving_target_id());
    }

    TEST_P(SBSearchDatabaseTest, AddGetObservatory)
    {
        const Observatory ztf{243.14022, 0.836325, +0.546877};
        const Observatory ldt{248.57749, 0.822887, 0.566916};
        const Observatory maunakea{204.5278, 0.94171, +0.33725};
        const Observatory paranal{289.59569, 0.909943, -0.414336};

        sbsdb->add_observatory("I41", ztf);
        sbsdb->add_observatory("G37", ldt);
        sbsdb->add_observatory("568", maunakea);
        sbsdb->add_observatory("309", paranal);

        Observatory obs = sbsdb->get_observatory("I41");
        EXPECT_EQ(obs, ztf);

        obs = sbsdb->get_observatory("G37");
        EXPECT_EQ(obs, ldt);

        obs = sbsdb->get_observatory("568");
        EXPECT_EQ(obs, maunakea);

        obs = sbsdb->get_observatory("309");
        EXPECT_EQ(obs, paranal);

        Observatories observatories = sbsdb->get_observatories();
        EXPECT_EQ(observatories["I41"], ztf);
        EXPECT_EQ(observatories["G37"], ldt);
        EXPECT_EQ(observatories["568"], maunakea);
        EXPECT_EQ(observatories["309"], paranal);

        sbsdb->remove_observatory("G37");
        EXPECT_THROW(sbsdb->get_observatory("G37"), ObservatoryError);

        EXPECT_THROW(sbsdb->add_observatory("I41", ztf), ObservatoryError);
    }

    TEST_P(SBSearchDatabaseTest, AddGetEphemeris)
    {
        MovingTarget encke{"2P"};
        Ephemeris eph{encke,
                      {{0, 10, 1, 0, 1, 0.1, 90, 0, 1, 180, 0, 0, 0, 10, -1},
                       {1, 11, 2, 0, 5, 0.5, 90, 1, 0, 0, 180, 30, 0, 20, 5},
                       {2, 12, 3, 0, 10, 1.0, 90, 2, 1, 90, 80, 90, 0, 30, 10}}};

        // The target is not in the database, so we expect an error
        EXPECT_THROW(sbsdb->add_ephemeris(eph), MovingTargetError);

        // Add the target, verify that the id was updated
        sbsdb->add_moving_target(encke);
        EXPECT_NE(encke.moving_target_id(), eph.target().moving_target_id());

        // Fix the target, and then we can add the ephemeris data
        eph.target(encke);
        sbsdb->add_ephemeris(eph);

        // Get the data back
        Ephemeris test;
        test = sbsdb->get_ephemeris(eph.target());
        EXPECT_EQ(test, eph);

        // Get a subset of data
        test = sbsdb->get_ephemeris(eph.target(), 0.5, 1.5);
        EXPECT_EQ(test, eph[1]);

        // This target does not match database copy
        MovingTarget wrong_id{"1P", eph.target().moving_target_id()};
        Ephemeris other{wrong_id, eph.data()};
        EXPECT_THROW(sbsdb->add_ephemeris(other), MovingTargetError);

        // Remove some data
        sbsdb->remove_ephemeris(eph.target(), 1.5, 10);
        test = sbsdb->get_ephemeris(eph.target());
        EXPECT_NE(test, eph);
        EXPECT_EQ(test, eph.slice(0, 2));

        // Remove all
        sbsdb->remove_ephemeris(eph.target());
        test = sbsdb->get_ephemeris(eph.target());
        EXPECT_EQ(test.num_vertices(), 0);
    }

    TEST_P(SBSearchDatabaseTest, AddGetObservation)
    {
        // nothing in the database yet
        EXPECT_EQ(sbsdb->count_observations(0, 10), 0);

        Observation obs("test source", "X05", "product", 0, 1, "0:0, 0:1, 1:1");
        // observation_id is not yet defined
        EXPECT_FALSE(obs.observation_id());

        // terms are not yet defined
        EXPECT_THROW(sbsdb->add_observations(obs), std::runtime_error);

        // update terms, add observation, now observation_id should be updated
        obs.terms(vector<string>{"asdf", "fdsa"});
        sbsdb->add_observations(obs);
        EXPECT_TRUE(obs.observation_id());
        EXPECT_EQ(sbsdb->count_observations(0, 10), 1);

        Observation retrieved = sbsdb->get_observation(obs.observation_id().value());
        EXPECT_TRUE(retrieved == obs);

        // edit the observation and update
        obs.terms(vector<string>{"a", "b", "c"});
        sbsdb->add_observations(obs);
        retrieved = sbsdb->get_observation(obs.observation_id().value());
        EXPECT_EQ(retrieved.terms(), vector<string>({"a", "b", "c"}));
        EXPECT_EQ(sbsdb->count_observations(0, 10), 1);

        // add a different source
        obs = Observation("another test source", "T05", "d", 4, 5, "0:0, 0:1, 1:1", "d e f");
        sbsdb->add_observations(obs);
        EXPECT_EQ(sbsdb->count_observations(0, 3), 1);
        EXPECT_EQ(sbsdb->count_observations(3, 10), 1);
        EXPECT_EQ(sbsdb->count_observations(0, 10), 2);
        EXPECT_EQ(sbsdb->count_observations("test source", 0, 10), 1);
        EXPECT_EQ(sbsdb->count_observations("another test source", 0, 10), 1);

        // try to get an observation that does not exist
        EXPECT_THROW(sbsdb->get_observation(-1), std::runtime_error);
    }

    TEST_P(SBSearchDatabaseTest, FindObservations)
    {
        Observation obs("test source", "X05", "a", 0, 1, "0:0, 0:1, 1:1", "a b c");
        sbsdb->add_observations(obs);

        obs = Observation("test source", "X05", "b", 1, 2, "0:0, 0:1, 1:1", "b c d");
        sbsdb->add_observations(obs);

        obs = Observation("test source", "X05", "c", 2, 3, "0:0, 0:1, 1:1", "c d e");
        sbsdb->add_observations(obs);

        obs = Observation("another test source", "T05", "d", 4, 5, "0:0, 0:1, 1:1", "d e f");
        sbsdb->add_observations(obs);

        // find observations matching term a
        set<int64_t> matches;
        matches = sbsdb->find_observation_ids(vector<string>{"a"});
        EXPECT_EQ(matches.size(), 1);

        // a or f
        matches = sbsdb->find_observation_ids(vector<string>{"a", "f"});
        EXPECT_EQ(matches.size(), 2);

        // c or f
        matches = sbsdb->find_observation_ids(vector<string>{"c", "f"});
        EXPECT_EQ(matches.size(), 4);

        // g
        matches = sbsdb->find_observation_ids(vector<string>{"g"});
        EXPECT_EQ(matches.size(), 0);

        // test observation time limits
        // start
        matches = sbsdb->find_observation_ids({"e"}, {.mjd_start = 2});
        EXPECT_EQ(matches.size(), 2);

        matches = sbsdb->find_observation_ids({"e"}, {.mjd_start = 3.5});
        EXPECT_EQ(matches.size(), 1);

        // stop
        matches = sbsdb->find_observation_ids({"e"}, {.mjd_stop = 1});
        EXPECT_EQ(matches.size(), 0);

        matches = sbsdb->find_observation_ids({"e"}, {.mjd_stop = 3});
        EXPECT_EQ(matches.size(), 1);

        matches = sbsdb->find_observation_ids({"e"}, {.mjd_stop = 5});
        EXPECT_EQ(matches.size(), 2);

        // start-stop
        matches = sbsdb->find_observation_ids({"e"}, {.mjd_start = 2, .mjd_stop = 2.5});
        EXPECT_EQ(matches.size(), 0);

        matches = sbsdb->find_observation_ids({"e"}, {.mjd_start = 2, .mjd_stop = 3});
        EXPECT_EQ(matches.size(), 1);

        matches = sbsdb->find_observation_ids({"e"}, {.mjd_start = 2.5, .mjd_stop = 4.5});
        EXPECT_EQ(matches.size(), 0);

        matches = sbsdb->find_observation_ids({"e"}, {.mjd_start = 3, .mjd_stop = 5});
        EXPECT_EQ(matches.size(), 1);

        // search by source
        matches = sbsdb->find_observation_ids({"b", "e"}, {.source = "test source"});
        EXPECT_EQ(matches.size(), 3);

        matches = sbsdb->find_observation_ids({"b", "e"}, {.source = "another test source"});
        EXPECT_EQ(matches.size(), 1);
    }

    TEST_P(SBSearchDatabaseTest, AddGetFound)
    {
        Observations observations({{"test source", "X05", "a", 0, 1, "0:0, 0:1, 1:1", "a b c"},
                                   {"test source", "X05", "b", 1, 2, "0:0, 0:1, 1:1", "b c d"}});
        sbsdb->add_observations(observations);

        MovingTarget encke("2P");
        sbsdb->add_moving_target(encke);

        Ephemeris eph{encke,
                      {{0, 10, 1, 0, 1, 0.1, 90, 0, 1, 180, 0, 0, 0, 10, -1},
                       {1, 11, 2, 0, 5, 0.5, 90, 1, 0, 0, 180, 30, 0, 20, 5},
                       {2, 12, 3, 0, 10, 1.0, 90, 2, 1, 90, 80, 90, 0, 30, 10}}};

        // these may not make sense, but it doesn't matter
        Founds founds;
        founds.append(Found(observations[0], eph[0]));
        founds.append(Found(observations[1], eph[1]));
        sbsdb->add_found(founds);

        founds = sbsdb->get_found(observations[0]);
        EXPECT_EQ(founds.size(), 1);
        EXPECT_EQ(founds[0], Found(observations[0], eph[0]));

        founds = sbsdb->get_found(observations[1]);
        EXPECT_EQ(founds.size(), 1);
        EXPECT_EQ(founds[0], Found(observations[1], eph[1]));

        founds = sbsdb->get_found(encke);
        EXPECT_EQ(founds.size(), 2);
        EXPECT_EQ(std::count(founds.begin(), founds.end(), Found(observations[0], eph[0])), 1);
        EXPECT_EQ(std::count(founds.begin(), founds.end(), Found(observations[1], eph[1])), 1);

        sbsdb->remove_found(founds);
        founds = sbsdb->get_found(observations[0]);
        EXPECT_EQ(founds.size(), 0);
        founds = sbsdb->get_found(observations[1]);
        EXPECT_EQ(founds.size(), 0);
        founds = sbsdb->get_found(encke);
        EXPECT_EQ(founds.size(), 0);
    }

    TEST_P(SBSearchDatabaseTest, ErrorIfClosed)
    {
        sbsdb->close();
        EXPECT_THROW(sbsdb->execute_sql("SELECT 1"), std::runtime_error);
    }

    INSTANTIATE_TEST_SUITE_P(Sqlite3AndPostgreSQL, SBSearchDatabaseTest, Values(&CreateSqlite3Database, &CreatePostgreSQLDatabase));
}