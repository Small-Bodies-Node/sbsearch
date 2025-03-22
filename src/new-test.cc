
#include <iostream>
#include <optional>
#include <vector>
#include <pqxx/pqxx>
#include <gtest/gtest.h>

#include "ephemeris.h"
#include "exceptions.h"
#include "indexer.h"
#include "moving_target.h"
#include "observation.h"
#include "observatory.h"
#include "sbsdb/add.h"
#include "sbsdb/get.h"
#include "sbsdb/remove.h"
#include "sbsdb/update.h"
#include "sbsdb/postgresql.h"

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

    TEST_F(SBSearchDatabaseTest, MovingTargetIO)
    {
        MovingTarget encke("2P");
        MovingTarget ceres("1");
        MovingTarget mercury("1", false);

        // add to the database, expect an updated moving_target_id
        EXPECT_FALSE(encke.moving_target_id());
        add::moving_target(db, encke);
        EXPECT_EQ(encke.moving_target_id(), 1);

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

        // add an alternate name and update encke
        encke.add_name("2P/Encke");
        update::moving_target(db, encke);
        test = get::moving_target(db, "2P");
        cerr << encke << " / " << test << endl;
        EXPECT_EQ(test.designation(), encke.designation());
        EXPECT_EQ(test.moving_target_id(), encke.moving_target_id());
        EXPECT_EQ(test.alternate_names(), encke.alternate_names());

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

    TEST_F(SBSearchDatabaseTest, AddGetObservatory)
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

    TEST_F(SBSearchDatabaseTest, AddGetObservations)
    {
        Observation obs("test source", "X05", "product", 0, 1, "0:0, 0:1, 1:1", "a b c");

        Observation result = add::observations(db, {obs})[0];
        EXPECT_FALSE(obs.observation_id());
        EXPECT_TRUE(result.observation_id());

        Observation retrieved = get::observation(db, result.observation_id().value());
        EXPECT_EQ(result, retrieved);
        EXPECT_EQ(result.terms(), retrieved.terms());
    }
}

// int main()
// {
//     Postgresql db("postgres:///lsst");

// auto count = db.get_one<optional<int>>("select count(*) from observations");
// cout << count.value_or(-1) << endl;
// auto fov = db.get_one<optional<string>>("select fov from observations limit 1");
// cout << fov.value_or("") << endl;

// auto observation = db.get_one<Observation>("select * from observations limit 1");
// cout << observation << endl;
// auto datum = db.get_one<Ephemeris::Datum>("select * from ephemerides limit 1");
// cout << Ephemeris(MovingTarget(), {datum}) << endl;

// auto observation_ids = db.get_many<int64_t>("select observation_id from observations limit 10");
// std::copy(observation_ids.begin(), observation_ids.end(), std::ostream_iterator<int64_t>(cout, "\n"));
// cout << endl;

// auto observations = Observations(db.get_many<Observation>("select * from observations limit 10"));
// cout << observations << endl;

// auto data = db.get_many<Ephemeris::Datum>("select * from ephemerides limit 10");
// cout << Ephemeris(MovingTarget(), data) << endl;

// auto target = get::moving_target(db, 1);
// cout << "Got by ID " << target << endl;

// target = get::moving_target(db, "Encke", true);
// cout << "Got by name " << target << endl;

// auto targets = get::all_moving_targets(db);
// cout << "All moving targets:\n  ";
// std::copy(targets.begin(), targets.end(), std::ostream_iterator<MovingTarget>(cout, "\n  "));
// cout << endl;

// observation = get::observation(db, observation_ids[5]);
// cout << observation << endl;

// observations = get::observations(db, observation_ids);
// cout << observations << endl;

// auto observatory = get::observatory(db, "X05");
// cout << observatory.longitude << endl;

// auto observatories = get::all_observatories(db);
// for (auto const &obs : observatories)
//     cout << obs.first << " " << obs.second.name << " " << obs.second.longitude << endl;

// auto sources = get::sources(db);
// std::copy(sources.begin(), sources.end(), std::ostream_iterator<string>(cout, "\n"));
// cout << endl;

// target = get::moving_target(db, 1);
// auto ephemeris = get::ephemeris(db, target, 0, 59590);
// cout << ephemeris << endl;

// observation = get::observation(db, 3);
// auto founds = get::found(db, observation);
// cout << founds << endl;

// observation = get::observation(db, 4);
// founds = get::found(db, observation);
// cout << founds << endl;

// target = get::moving_target(db, 1);
// founds = get::found(db, target);
// cout << founds << endl;

// target = get::moving_target(db, 5);
// founds = get::found(db, target);
// cout << founds << endl;

// target = MovingTarget("11P");
// target.add_name("TSL");
// add::target()
// }
