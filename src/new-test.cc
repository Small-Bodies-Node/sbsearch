
#include <iostream>
#include <optional>
#include <unordered_set>
#include <vector>
#include <pqxx/pqxx>
#include <gtest/gtest.h>

#include "ephemeris.h"
#include "exceptions.h"
#include "indexer.h"
#include "moving_target.h"
#include "observation.h"
#include "observatory.h"
#include "sbsdb/sbsdb.h"

using namespace sbsearch;
using namespace sbsearch::sbsdb;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace sbsearch::testing
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

    TEST_F(SBSearchDatabaseTest, SetupTables)
    {
        EXPECT_EQ(db.get_one<int>("SELECT COUNT(*) FROM configuration"), 5);
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

    TEST_F(SBSearchDatabaseTest, EphemerisIO)
    {
        MovingTarget encke{"2P"};
        Ephemeris eph{encke,
                      {{0, 10, 1, 0, 1, 0.1, 90, 0, 1, 180, 0, 0, 0, 10, -1},
                       {1, 11, 2, 0, 5, 0.5, 90, 1, 0, 0, 180, 30, 0, 20, 5},
                       {2, 12, 3, 0, 10, 1.0, 90, 2, 1, 90, 80, 90, 0, 30, 10}}};

        // The target is not in the database, so we expect an error
        EXPECT_THROW(add::ephemeris(db, eph), MovingTargetError);

        // Add the target, verify that the id was updated
        add::moving_target(db, encke);
        EXPECT_NE(encke.moving_target_id(), eph.target().moving_target_id());

        // Fix the target, and then we can add the ephemeris data
        eph.target(encke);
        add::ephemeris(db, eph);

        // Get the data back
        Ephemeris test;
        test = get::ephemeris(db, eph.target());
        EXPECT_EQ(test, eph);

        // Get a subset of data
        test = get::ephemeris(db, eph.target(), 0.5, 1.5);
        EXPECT_EQ(test, eph[1]);

        // This target does not match database copy
        MovingTarget wrong_id{"1P", eph.target().moving_target_id()};
        Ephemeris other{wrong_id, eph.data()};
        EXPECT_THROW(add::ephemeris(db, other), MovingTargetError);

        // Remove some data
        remove::ephemeris(db, eph.target(), 1.5, 10);
        test = get::ephemeris(db, eph.target());
        EXPECT_NE(test, eph);
        EXPECT_EQ(test, eph.slice(0, 2));

        // Remove all
        remove::ephemeris(db, eph.target());
        test = get::ephemeris(db, eph.target());
        EXPECT_EQ(test.num_vertices(), 0);
    }

    TEST_F(SBSearchDatabaseTest, EphemerisDateRange)
    {
        MovingTarget encke{"2P"};
        EXPECT_THROW(get::ephemeris_date_range(db, encke), MovingTargetError); // need obsid

        add::moving_target(db, encke);
        auto range = get::ephemeris_date_range(db, encke);
        EXPECT_FALSE(range.first); // no ephemeris data
        EXPECT_FALSE(range.second);

        Ephemeris encke_eph{encke,
                            {{0, 10, 1, 0, 1, 0.1, 90, 0, 1, 180, 0, 0, 0, 10, -1},
                             {1, 11, 2, 0, 5, 0.5, 90, 1, 0, 0, 180, 30, 0, 20, 5},
                             {2, 12, 3, 0, 10, 1.0, 90, 2, 1, 90, 80, 90, 0, 30, 10}}};
        add::ephemeris(db, encke_eph);

        MovingTarget tempel1("9P");
        add::moving_target(db, tempel1);
        Ephemeris tempel1_eph{tempel1,
                              {{1.5, 10, 1, 0, 1, 0.1, 90, 0, 1, 180, 0, 0, 0, 10, -1},
                               {2.5, 11, 2, 0, 5, 0.5, 90, 1, 0, 0, 180, 30, 0, 20, 5},
                               {3.5, 12, 3, 0, 10, 1.0, 90, 2, 1, 90, 80, 90, 0, 30, 10}}};
        add::ephemeris(db, tempel1_eph);

        range = get::ephemeris_date_range(db, encke);
        EXPECT_EQ(range.first, 0);
        EXPECT_EQ(range.second, 2);

        range = get::ephemeris_date_range(db, tempel1);
        EXPECT_EQ(range.first, 1.5);
        EXPECT_EQ(range.second, 3.5);
    }

    TEST_F(SBSearchDatabaseTest, FoundIO)
    {
        Observations observations({{"test source", "X05", "a", 0, 1, "0:0, 0:1, 1:1", "a b c"},
                                   {"test source", "X05", "b", 1, 2, "0:0, 0:1, 1:1", "b c d"}});
        add::observations(db, observations);

        MovingTarget encke("2P");
        add::moving_target(db, encke);

        Ephemeris eph{encke,
                      {{0, 10, 1, 0, 1, 0.1, 90, 0, 1, 180, 0, 0, 0, 10, -1},
                       {1, 11, 2, 0, 5, 0.5, 90, 1, 0, 0, 180, 30, 0, 20, 5},
                       {2, 12, 3, 0, 10, 1.0, 90, 2, 1, 90, 80, 90, 0, 30, 10}}};

        // these may not make sense, but it doesn't matter
        Founds founds;
        founds.append(Found(observations[0], eph.segment(0)));
        founds.append(Found(observations[1], eph.segment(1)));
        add::found(db, founds);

        founds = get::found(db, observations[0]);
        EXPECT_EQ(founds.size(), 1);
        EXPECT_EQ(founds[0], Found(observations[0], eph[0]));

        founds = get::found(db, observations[1]);
        EXPECT_EQ(founds.size(), 1);
        EXPECT_EQ(founds[0], Found(observations[1], eph[1]));

        founds = get::found(db, encke);
        EXPECT_EQ(founds.size(), 2);
        EXPECT_EQ(std::count(founds.begin(), founds.end(), Found(observations[0], eph[0])), 1);
        EXPECT_EQ(std::count(founds.begin(), founds.end(), Found(observations[1], eph[1])), 1);

        remove::found(db, founds);
        founds = get::found(db, observations[0]);
        EXPECT_EQ(founds.size(), 0);
        founds = get::found(db, observations[1]);
        EXPECT_EQ(founds.size(), 0);
        founds = get::found(db, encke);
        EXPECT_EQ(founds.size(), 0);
    }

    TEST_F(SBSearchDatabaseTest, IndexerOptionsIO)
    {
        auto options = get::indexer_options(db);

        // expect the defaults
        EXPECT_EQ(options.max_spatial_index_cells(), 8);
        EXPECT_EQ(options.max_spatial_level(), 12);
        EXPECT_EQ(options.min_spatial_level(), 4);
        EXPECT_EQ(options.temporal_resolution(), 1);

        // change everything and verify
        Indexer::Options new_options;
        new_options.max_spatial_index_cells(10);
        new_options.max_spatial_level(11);
        new_options.min_spatial_level(6);
        new_options.temporal_resolution(2);
        update::indexer_options(db, new_options);
        EXPECT_TRUE(new_options == get::indexer_options(db));
    }

    TEST_F(SBSearchDatabaseTest, MovingTargetIO)
    {
        MovingTarget encke("2P");
        MovingTarget ceres("1");
        MovingTarget mercury("1", false);

        // verify encke should fail: not yet in the database
        EXPECT_THROW(verify::moving_target(db, encke), MovingTargetError);

        // add to the database, expect an updated moving_target_id
        EXPECT_FALSE(encke.moving_target_id());
        add::moving_target(db, encke);
        EXPECT_EQ(encke.moving_target_id(), 1);

        // verify encke should now pass
        EXPECT_NO_THROW(verify::moving_target(db, encke));

        // add another, it should be 2
        add::moving_target(db, ceres);
        EXPECT_EQ(ceres.moving_target_id(), 2);

        // and this should be 3
        add::moving_target(db, mercury);
        EXPECT_EQ(mercury.moving_target_id(), 3);

        // get them from the database
        MovingTarget test;
        test = get::moving_target(db, 1);
        EXPECT_EQ(test.designation(), encke.designation());
        EXPECT_EQ(test.moving_target_id(), encke.moving_target_id());
        EXPECT_EQ(test.alternate_names(), encke.alternate_names());
        EXPECT_EQ(test.small_body(), encke.small_body());

        test = get::moving_target(db, "2P");
        EXPECT_EQ(test.designation(), encke.designation());
        EXPECT_EQ(test.moving_target_id(), encke.moving_target_id());
        EXPECT_EQ(test.alternate_names(), encke.alternate_names());
        EXPECT_EQ(test.small_body(), encke.small_body());

        test = get::moving_target(db, "1");
        EXPECT_EQ(test.designation(), ceres.designation());
        EXPECT_EQ(test.moving_target_id(), ceres.moving_target_id());
        EXPECT_EQ(test.alternate_names(), ceres.alternate_names());
        EXPECT_EQ(test.small_body(), ceres.small_body());

        test = get::moving_target(db, 2);
        EXPECT_EQ(test.designation(), ceres.designation());
        EXPECT_EQ(test.moving_target_id(), ceres.moving_target_id());
        EXPECT_EQ(test.alternate_names(), ceres.alternate_names());
        EXPECT_EQ(test.small_body(), ceres.small_body());

        test = get::moving_target(db, "1", false);
        EXPECT_EQ(test.designation(), mercury.designation());
        EXPECT_EQ(test.moving_target_id(), mercury.moving_target_id());
        EXPECT_EQ(test.alternate_names(), mercury.alternate_names());
        EXPECT_EQ(test.small_body(), mercury.small_body());

        test = get::moving_target(db, 3);
        EXPECT_EQ(test.designation(), mercury.designation());
        EXPECT_EQ(test.moving_target_id(), mercury.moving_target_id());
        EXPECT_EQ(test.alternate_names(), mercury.alternate_names());
        EXPECT_EQ(test.small_body(), mercury.small_body());

        // add an alternate name to encke and it no longer matches the db copy
        encke.add_name("2P/Encke");
        EXPECT_THROW(verify::moving_target(db, encke), MovingTargetError);

        // update and it should pass
        update::moving_target(db, encke);
        EXPECT_NO_THROW(verify::moving_target(db, encke));

        // try getting Encke via alt name
        test = get::moving_target(db, "2P/Encke");
        EXPECT_EQ(test.designation(), encke.designation());
        EXPECT_EQ(test.moving_target_id(), encke.moving_target_id());
        EXPECT_EQ(test.alternate_names(), encke.alternate_names());

        // add a few names to Ceres
        vector<string> names{"(1) Ceres", "Ceres", "A801 AA"};
        ceres.add_names(names.begin(), names.end());
        update::moving_target(db, ceres);
        test = get::moving_target(db, "A801 AA");
        EXPECT_EQ(test.designation(), ceres.designation());
        EXPECT_EQ(test.moving_target_id(), ceres.moving_target_id());
        EXPECT_EQ(test.alternate_names(), ceres.alternate_names());
        EXPECT_EQ(test.alternate_names().size(), 3);

        // and Mercury
        names = vector<string>{"Hermes", "Nabu"};
        mercury.add_names(names.begin(), names.end());
        update::moving_target(db, mercury);
        test = get::moving_target(db, "Nabu", false);
        EXPECT_EQ(test.designation(), mercury.designation());
        EXPECT_EQ(test.moving_target_id(), mercury.moving_target_id());
        EXPECT_EQ(test.alternate_names(), mercury.alternate_names());
        EXPECT_EQ(test.alternate_names().size(), 2);

        // Try to add an object that already exists.
        MovingTarget halley("1P");
        halley.moving_target_id(1);
        EXPECT_THROW(add::moving_target(db, halley), MovingTargetError);

        MovingTarget duplicate_ceres("1");
        EXPECT_THROW(add::moving_target(db, duplicate_ceres), MovingTargetError);

        // Removing an object that does not exist is just a warning.
        MovingTarget new_comet("1000P");
        new_comet.moving_target_id(9123);
        remove::moving_target(db, new_comet);

        // Get objects that do not exist
        EXPECT_THROW(get::moving_target(db, 1000), MovingTargetError);
        EXPECT_EQ(get::moving_target(db, "asdf"), MovingTarget("asdf"));
        EXPECT_FALSE(get::moving_target(db, "Nabu", true).moving_target_id());
        EXPECT_FALSE(get::moving_target(db, "Ceres", false).moving_target_id());
    }

    TEST_F(SBSearchDatabaseTest, ObservationsIO)
    {
        // nothing in the database yet
        EXPECT_EQ(count::observations(db, 0, 10), 0);

        Observations obs{{"test source", "X05", "product", 0, 1, "0:0, 0:1, 1:1"}};

        // verify that observation_id is not yet defined
        EXPECT_NO_THROW(verify::observations(obs, false, false));

        // but required terms are not yet defined
        EXPECT_THROW(verify::observations(obs, false, true), ObservationError);

        // update terms, add observation, and check updated observation_id
        obs[0].terms(vector<string>{"asdf", "fdsa"});
        add::observations(db, obs);
        EXPECT_TRUE(obs[0].observation_id());
        EXPECT_EQ(count::observations(db, 0, 10), 1);

        // get that observation from the database and check that it matches
        Observations retrieved = get::observations(db, obs.observation_ids());
        EXPECT_TRUE(retrieved[0] == obs[0]);

        // before updating, verify that observation_ids and terms are defined
        EXPECT_NO_THROW(verify::observations(obs, true, true));

        // edit the observation and update
        obs[0].terms(vector<string>{"a", "b", "c"});
        update::observations(db, {obs});
        retrieved = get::observations(db, obs.observation_ids());
        EXPECT_EQ(retrieved[0].terms(), vector<string>({"a", "b", "c"}));
        EXPECT_EQ(count::observations(db, 0, 10), 1);

        // add a different source
        obs = Observations{{"another test source", "T05", "d", 4, 5, "0:0, 0:1, 1:1", "d e f"}};
        add::observations(db, obs);
        EXPECT_EQ(count::observations(db, 0, 3), 1);
        EXPECT_EQ(count::observations(db, 3, 10), 1);
        EXPECT_EQ(count::observations(db, 0, 10), 2);
        EXPECT_EQ(count::observations(db, "test source", 0, 10), 1);
        EXPECT_EQ(count::observations(db, "another test source", 0, 10), 1);

        // remove an observation and try to get it from the database.
        auto ids = obs.observation_ids();
        remove::observations(db, obs);
        EXPECT_FALSE(obs[0].observation_id());
        EXPECT_THROW(get::observations(db, ids), std::runtime_error);
    }

    TEST_F(SBSearchDatabaseTest, AllObservationsFOV)
    {
        Observations observations({
            Observation("test source 1", "X05", "product1", 0, 1, "0:0, 0:1, 1:1", "a b c"),
            Observation("test source 2", "568", "product2", 1, 2, "0:1, 0:2, 1:2", "b c d"),
            Observation("test source 1", "X05", "product3", 2, 3, "0:2, 0:3, 1:3", "c d e"),
            Observation("test source 2", "568", "product4", 3, 4, "0:3, 0:4, 1:4", "d e f"),
        });
        add::observations(db, observations);

        auto rows = get::all_observations_fov(db, 2, 0);
        EXPECT_EQ(rows.size(), 2);
        EXPECT_EQ(rows[0].first, observations[0].observation_id());
        EXPECT_EQ(rows[0].second, "0:0, 0:1, 1:1");
        EXPECT_EQ(rows[1].first, observations[1].observation_id());
        EXPECT_EQ(rows[1].second, "0:1, 0:2, 1:2");

        rows = get::all_observations_fov(db, 2, 2);
        EXPECT_EQ(rows.size(), 2);
        EXPECT_EQ(rows[0].first, observations[2].observation_id());
        EXPECT_EQ(rows[0].second, "0:2, 0:3, 1:3");
        EXPECT_EQ(rows[1].first, observations[3].observation_id());
        EXPECT_EQ(rows[1].second, "0:3, 0:4, 1:4");

        rows = get::all_observations_fov(db, 4, 1);
        EXPECT_EQ(rows.size(), 3);
        EXPECT_EQ(rows[0].first, observations[1].observation_id());
        EXPECT_EQ(rows[0].second, "0:1, 0:2, 1:2");
        EXPECT_EQ(rows[1].first, observations[2].observation_id());
        EXPECT_EQ(rows[1].second, "0:2, 0:3, 1:3");
        EXPECT_EQ(rows[2].first, observations[3].observation_id());
        EXPECT_EQ(rows[2].second, "0:3, 0:4, 1:4");

        rows = get::all_observations_fov(db, 4, 4);
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
            observations[i].terms("asdf fdsa");

        add::observations(db, observations);

        auto drange = get::observations_date_range(db);
        EXPECT_EQ(drange.first, 0);
        EXPECT_EQ(drange.second, 4);

        drange = get::observations_date_range(db, "test source 1");
        EXPECT_EQ(drange.first, 0);
        EXPECT_EQ(drange.second, 3);

        drange = get::observations_date_range(db, "test source 2");
        EXPECT_EQ(drange.first.value(), 1);
        EXPECT_EQ(drange.second.value(), 4);

        // null for no observations
        drange = get::observations_date_range(db, "test source 3");
        EXPECT_FALSE(drange.first);
        EXPECT_FALSE(drange.second);
    }

    TEST_F(SBSearchDatabaseTest, ObservatoryIO)
    {
        const Observatory ztf{243.14022, 0.836325, +0.546877, "I41"};
        const Observatory ldt{248.57749, 0.822887, 0.566916, "G37"};
        const Observatory maunakea{204.5278, 0.94171, +0.33725, "568"};
        const Observatory paranal{289.59569, 0.909943, -0.414336, "309"};

        add::observatory(db, ztf);
        add::observatory(db, ldt);
        add::observatory(db, maunakea);
        add::observatory(db, paranal);

        Observatory obs = get::observatory(db, "I41");
        EXPECT_EQ(obs, ztf);

        obs = get::observatory(db, "G37");
        EXPECT_EQ(obs, ldt);

        obs = get::observatory(db, "568");
        EXPECT_EQ(obs, maunakea);

        obs = get::observatory(db, "309");
        EXPECT_EQ(obs, paranal);

        Observatories observatories = get::all_observatories(db);
        EXPECT_EQ(observatories["I41"], ztf);
        EXPECT_EQ(observatories["G37"], ldt);
        EXPECT_EQ(observatories["568"], maunakea);
        EXPECT_EQ(observatories["309"], paranal);

        remove::observatory(db, "G37");
        EXPECT_THROW(get::observatory(db, "G37"), ObservatoryError);

        EXPECT_THROW(add::observatory(db, ztf), ObservatoryError);
    }

    TEST_F(SBSearchDatabaseTest, FindObservations)
    {
        Observations obs({{"test source", "X05", "a", 0, 1, "0:0, 0:1, 1:1", "a b c"},
                          {"test source", "X05", "b", 1, 2, "0:0, 0:1, 1:1", "b c d"},
                          {"test source", "X05", "c", 2, 3, "0:0, 0:1, 1:1", "c d e"},
                          {"another test source", "T05", "d", 4, 5, "0:0, 0:1, 1:1", "d e f"}});
        add::observations(db, obs);

        // find observations matching term a
        std::unordered_set<int64_t> matches;
        matches = find::observations(db, vector<string>{"a"});
        EXPECT_EQ(matches.size(), 1);

        // a or f
        matches = find::observations(db, vector<string>{"a", "f"});
        EXPECT_EQ(matches.size(), 2);

        // c or f
        matches = find::observations(db, vector<string>{"c", "f"});
        EXPECT_EQ(matches.size(), 4);

        // g
        matches = find::observations(db, vector<string>{"g"});
        EXPECT_EQ(matches.size(), 0);

        // test observation time limits
        // start
        matches = find::observations(db, {"e"}, {.mjd_start = 2});
        EXPECT_EQ(matches.size(), 2);

        matches = find::observations(db, {"e"}, {.mjd_start = 3.5});
        EXPECT_EQ(matches.size(), 1);

        // stop
        matches = find::observations(db, {"e"}, {.mjd_stop = 1});
        EXPECT_EQ(matches.size(), 0);

        matches = find::observations(db, {"e"}, {.mjd_stop = 3});
        EXPECT_EQ(matches.size(), 1);

        matches = find::observations(db, {"e"}, {.mjd_stop = 5});
        EXPECT_EQ(matches.size(), 2);

        // start-stop
        matches = find::observations(db, {"e"}, {.mjd_start = 2, .mjd_stop = 2.5});
        EXPECT_EQ(matches.size(), 0);

        matches = find::observations(db, {"e"}, {.mjd_start = 2, .mjd_stop = 3});
        EXPECT_EQ(matches.size(), 1);

        matches = find::observations(db, {"e"}, {.mjd_start = 2.5, .mjd_stop = 4.5});
        EXPECT_EQ(matches.size(), 0);

        matches = find::observations(db, {"e"}, {.mjd_start = 3, .mjd_stop = 5});
        EXPECT_EQ(matches.size(), 1);

        // search by source
        matches = find::observations(db, {"b", "e"}, {.source = "test source"});
        EXPECT_EQ(matches.size(), 3);

        matches = find::observations(db, {"b", "e"}, {.source = "another test source"});
        EXPECT_EQ(matches.size(), 1);
    }

}
