#include "config.h"

#include <iterator>
#include <string>
#include <vector>
#include <boost/json.hpp>
#include <gtest/gtest.h>

#include "date.h"
#include "found.h"
#include "json.h"
#include "observation.h"
#include "ephemeris/ephemeris.h"
#include "ephemeris/interpolate.h"

using sbsearch::Observation;
using sbsearch::ephemeris::Ephemeris;
using std::string;
using std::vector;

namespace sbsearch::testing
{
    TEST(FoundTests, StreamInsertOperator)
    {
        MovingTarget encke{"2P"};

        Observations observations(
            {{"test source", "I41", "a", 59252.01, 59252.019, "1:3, 2:3, 2:4, 1:4", {}, 1},
             {"test source", "I41", "b", 59252.02, 59252.029, "2:3, 3:3, 3:4, 2:4", {}, 2}});
        observations[0].meta("test=0.9");

        Ephemeris eph(encke, {{59252.01, 10.01, 0.0, 3.5, 375, 90, 0, 0, 0, 1, 1, 0},
                              {59252.02, 10.02, 1.5, 3.5, 375, 90, 0, 0, 0, 1, 1, 0},
                              {59252.03, 10.03, 2.5, 3.5, 375, 90, 0, 0, 0, 1, 1, 0},
                              {59252.04, 10.04, 3.5, 3.5, 375, 90, 0, 0, 0, 1, 1, 0}});
        Founds founds;

        founds.data.emplace_back(observations[0], eph.segment(0));
        founds.data.emplace_back(observations[1], ephemeris::interpolate(eph, observations[1].mjd_mid()));

        // Should be two found observations
        EXPECT_EQ(founds.size(), 2);

        std::stringstream stream;
        stream.str("");
        stream << founds;
        EXPECT_EQ(
            stream.str(),
            "observation_id       source  product_id  observatory     mjd_start      mjd_stop  exposure  moving_target_id  designation  small_body           mjd       tmtp        ra       dec      mu  mu_theta      rh   delta  phase  selong  true_anomaly  sangle  vangle  unc_a  unc_b  unc_th  vmag\n"
            "--------------  -----------  ----------  -----------  ------------  ------------  --------  ----------------  -----------  ----------  ------------  ---------  --------  --------  ------  --------  ------  ------  -----  ------  ------------  ------  ------  -----  -----  ------  ----\n"
            "             1  test source           a          I41  59252.010000  59252.019000   777.600              null           2P        true  59252.014500  10.014500  0.674996  3.500296  375.00    90.000  1.0000  1.0000  0.000    null          null    null    null  0.000  0.000   0.000  null\n"
            "             2  test source           b          I41  59252.020000  59252.029000   777.600              null           2P        true  59252.024500  10.024500  1.981954  3.499942  375.00    90.000  1.0000  1.0000  0.000    null          null    null    null  0.000  0.000   0.000  null\n");

        stream.str("");
        founds.observation_format.show_fov = true;
        stream << founds;
        EXPECT_EQ(
            stream.str(),
            "observation_id       source  product_id  observatory     mjd_start      mjd_stop  exposure                 fov  moving_target_id  designation  small_body           mjd       tmtp        ra       dec      mu  mu_theta      rh   delta  phase  selong  true_anomaly  sangle  vangle  unc_a  unc_b  unc_th  vmag\n"
            "--------------  -----------  ----------  -----------  ------------  ------------  --------  ------------------  ----------------  -----------  ----------  ------------  ---------  --------  --------  ------  --------  ------  ------  -----  ------  ------------  ------  ------  -----  -----  ------  ----\n"
            "             1  test source           a          I41  59252.010000  59252.019000   777.600  1:3, 2:3, 2:4, 1:4              null           2P        true  59252.014500  10.014500  0.674996  3.500296  375.00    90.000  1.0000  1.0000  0.000    null          null    null    null  0.000  0.000   0.000  null\n"
            "             2  test source           b          I41  59252.020000  59252.029000   777.600  2:3, 3:3, 3:4, 2:4              null           2P        true  59252.024500  10.024500  1.981954  3.499942  375.00    90.000  1.0000  1.0000  0.000    null          null    null    null  0.000  0.000   0.000  null\n");

        stream.str("");
        founds.observation_format.show_fov = false;
        founds.ephemeris_format.date = Date::Format::Calendar;
        stream << founds;
        EXPECT_EQ(
            stream.str(),
            "observation_id       source  product_id  observatory     mjd_start      mjd_stop  exposure  moving_target_id  designation  small_body                 date       tmtp        ra       dec      mu  mu_theta      rh   delta  phase  selong  true_anomaly  sangle  vangle  unc_a  unc_b  unc_th  vmag\n"
            "--------------  -----------  ----------  -----------  ------------  ------------  --------  ----------------  -----------  ----------  -------------------  ---------  --------  --------  ------  --------  ------  ------  -----  ------  ------------  ------  ------  -----  -----  ------  ----\n"
            "             1  test source           a          I41  59252.010000  59252.019000   777.600              null           2P        true  2021-02-07 00:20:53  10.014500  0.674996  3.500296  375.00    90.000  1.0000  1.0000  0.000    null          null    null    null  0.000  0.000   0.000  null\n"
            "             2  test source           b          I41  59252.020000  59252.029000   777.600              null           2P        true  2021-02-07 00:35:17  10.024500  1.981954  3.499942  375.00    90.000  1.0000  1.0000  0.000    null          null    null    null  0.000  0.000   0.000  null\n");
    }

    TEST(FoundTests, FoundAsJSON)
    {
        MovingTarget encke{"2P"};

        Observation obs("test source", "I41", "a", 59252.01, 59252.019, "1:3, 2:3, 2:4, 1:4", {}, 1);
        obs.meta("poor seeing");
        obs.format.show_fov = true;
        obs.format.show_meta = true;

        Ephemeris eph(encke, {{59252.01, 10.01, 0, 3.5, 375, 90, 0, 0, 0, 1, 1, 0},
                              {59252.02, 10.02, 1.5, 3.5, 375, 90, 0, 0, 0, 1, 1, 0}});
        Found found(obs, eph);

        json::object obj = json::value_from(found).as_object();

        EXPECT_EQ(obj["source"], "test source");
        EXPECT_EQ(obj["observatory"], "I41");
        EXPECT_EQ(obj["product_id"], "a");
        EXPECT_EQ(obj["mjd_start"], 59252.01);
        EXPECT_EQ(obj["mjd_stop"], 59252.019);
        EXPECT_EQ(obj["fov"], "1:3, 2:3, 2:4, 1:4");
        EXPECT_EQ(obj["meta"], "poor seeing");
        EXPECT_FLOAT_EQ(*obj["mjd"].if_double(), 59252.014500);
        EXPECT_FLOAT_EQ(*obj["tmtp"].if_double(), 10.0145);
        EXPECT_FLOAT_EQ(*obj["ra"].if_double(), 0.67499578);
        EXPECT_FLOAT_EQ(*obj["dec"].if_double(), 3.5002961);
        EXPECT_FLOAT_EQ(*obj["mu"].if_double(), 375);
        EXPECT_FLOAT_EQ(*obj["mu_theta"].if_double(), 90);
        EXPECT_EQ(obj["unc_a"], 0.);
        EXPECT_EQ(obj["unc_b"], 0.);
        EXPECT_EQ(obj["unc_theta"], 0.);
        EXPECT_EQ(obj["rh"], 1.);
        EXPECT_EQ(obj["delta"], 1.);
        EXPECT_EQ(obj["phase"], 0.);
        EXPECT_EQ(obj["selong"], nullptr);
        EXPECT_EQ(obj["true_anomaly"], nullptr);
        EXPECT_EQ(obj["sangle"], nullptr);
        EXPECT_EQ(obj["vangle"], nullptr);
        EXPECT_EQ(obj["vmag"], nullptr);

        found.ephemeris.format.date = Date::Format::Calendar;
        obj = json::value_from(found).as_object();
        EXPECT_EQ(obj["date"], "2021-02-07 00:20:53");
    }

    TEST(FoundsTest, FoundsAsJSON)
    {
        MovingTarget encke{"2P"};

        Observations observations({Observation("test source", "I41", "a", 59252.01, 59252.019, "1:3, 2:3, 2:4, 1:4", {}, 1, {}, "meta a_1"),
                                   Observation("test source", "I41", "b", 59252.02, 59252.029, "2:3, 3:3, 3:4, 2:4", {}, 2, {}, "meta b_2")});

        Ephemeris eph(encke, {{59252.01, 10.01, 0.0, 3.5, 375, 90, 0, 0, 0, 1, 1, 0},
                              {59252.02, 10.02, 1.5, 3.5, 375, 90, 0, 0, 0, 1, 1, 0},
                              {59252.03, 10.03, 2.5, 3.5, 375, 90, 0, 0, 0, 1, 1, 0},
                              {59252.04, 10.04, 3.5, 3.5, 375, 90, 0, 0, 0, 1, 1, 0}});
        Founds founds;

        founds.data.emplace_back(observations[0], eph.segment(0));
        founds.append(Found(observations[1], ephemeris::interpolate(eph, observations[1].mjd_mid())));
        founds.observation_format.show_fov = false;
        founds.observation_format.show_meta = false;
        founds.ephemeris_format.date = Date::Format::MJD;

        json::array array = json::value_from(founds).as_array();

        json::object obj = *array.at(0).if_object();
        EXPECT_EQ(obj["source"], "test source");
        EXPECT_EQ(obj["observatory"], "I41");
        EXPECT_EQ(obj["product_id"], "a");
        EXPECT_EQ(obj["mjd_start"], 59252.01);
        EXPECT_EQ(obj["mjd_stop"], 59252.019);
        EXPECT_TRUE(obj["fov"].is_null());
        EXPECT_TRUE(obj["meta"].is_null());
        EXPECT_FLOAT_EQ(*obj["mjd"].if_double(), 59252.014500);
        EXPECT_FLOAT_EQ(*obj["tmtp"].if_double(), 10.0145);
        EXPECT_FLOAT_EQ(*obj["ra"].if_double(), 0.67499578);
        EXPECT_FLOAT_EQ(*obj["dec"].if_double(), 3.500296);
        EXPECT_FLOAT_EQ(*obj["mu"].if_double(), 375);
        EXPECT_FLOAT_EQ(*obj["mu_theta"].if_double(), 90);
        EXPECT_EQ(obj["unc_a"], 0.);
        EXPECT_EQ(obj["unc_b"], 0.);
        EXPECT_EQ(obj["unc_theta"], 0.);
        EXPECT_EQ(obj["rh"], 1.);
        EXPECT_EQ(obj["delta"], 1.);
        EXPECT_EQ(obj["phase"], 0.);
        EXPECT_EQ(obj["selong"], nullptr);
        EXPECT_EQ(obj["true_anomaly"], nullptr);
        EXPECT_EQ(obj["sangle"], nullptr);
        EXPECT_EQ(obj["vangle"], nullptr);
        EXPECT_EQ(obj["vmag"], nullptr);

        obj = *array.at(1).if_object();
        EXPECT_EQ(obj["source"], "test source");
        EXPECT_EQ(obj["observatory"], "I41");
        EXPECT_EQ(obj["product_id"], "b");
        EXPECT_EQ(obj["mjd_start"], 59252.02);
        EXPECT_EQ(obj["mjd_stop"], 59252.029);
        EXPECT_TRUE(obj["fov"].is_null());
        EXPECT_TRUE(obj["meta"].is_null());
        EXPECT_FLOAT_EQ(*obj["mjd"].if_double(), 59252.024500);
        EXPECT_FLOAT_EQ(*obj["tmtp"].if_double(), 10.0245);
        EXPECT_FLOAT_EQ(*obj["ra"].if_double(), 1.9819541);
        EXPECT_FLOAT_EQ(*obj["dec"].if_double(), 3.4999416);
        EXPECT_FLOAT_EQ(*obj["mu"].if_double(), 375);
        EXPECT_FLOAT_EQ(*obj["mu_theta"].if_double(), 90);
        EXPECT_EQ(obj["unc_a"], 0.);
        EXPECT_EQ(obj["unc_b"], 0.);
        EXPECT_EQ(obj["unc_theta"], 0.);
        EXPECT_EQ(obj["rh"], 1.);
        EXPECT_EQ(obj["delta"], 1.);
        EXPECT_EQ(obj["phase"], 0.);
        EXPECT_EQ(obj["selong"], nullptr);
        EXPECT_EQ(obj["true_anomaly"], nullptr);
        EXPECT_EQ(obj["sangle"], nullptr);
        EXPECT_EQ(obj["vangle"], nullptr);
        EXPECT_EQ(obj["vmag"], nullptr);

        // test formatting
        founds.observation_format.show_fov = true;
        founds.observation_format.show_meta = true;
        founds.observation_format.date = Date::Format::Calendar;
        founds.ephemeris_format.date = Date::Format::Calendar;

        array = json::value_from(founds).as_array();
        obj = *array.at(0).if_object();
        EXPECT_EQ(obj["mjd_start"], "2021-02-07 00:14:24");
        EXPECT_EQ(obj["mjd_stop"], "2021-02-07 00:27:22");
        EXPECT_EQ(obj["date"], "2021-02-07 00:20:53");
        EXPECT_EQ(obj["fov"], "1:3, 2:3, 2:4, 1:4");
        EXPECT_EQ(obj["meta"], "meta a_1");

        obj = *array.at(1).if_object();
        EXPECT_EQ(obj["mjd_start"], "2021-02-07 00:28:48");
        EXPECT_EQ(obj["mjd_stop"], "2021-02-07 00:41:46");
        EXPECT_EQ(obj["date"], "2021-02-07 00:35:17");
        EXPECT_EQ(obj["fov"], "2:3, 3:3, 3:4, 2:4");
        EXPECT_EQ(obj["meta"], "meta b_2");
    }
}