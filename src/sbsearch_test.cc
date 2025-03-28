#include "config.h"

// #include <cstdio>
// #include <sstream>
#include <string>
#include <vector>
// #include <unistd.h>
#include <gtest/gtest.h>
#include <s2/s2latlng.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>

#include "ephemeris.h"
#include "exceptions.h"
#include "indexer.h"
#include "observation.h"
#include "sbsdb/sbsdb.h"
#include "sbsearch.h"

using namespace sbsearch;
using std::vector;

class SBSearchTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        sbs.db()->execute(R"(
            DROP TABLE IF EXISTS observation_terms_index;
            DROP TABLE IF EXISTS configuration;
            DROP TABLE IF EXISTS found;
            DROP TABLE IF EXISTS ephemerides;
            DROP TABLE IF EXISTS observatories;
            DROP TABLE IF EXISTS moving_targets;
            DROP TABLE IF EXISTS observations;
        )");
        sbs.db()->setup_tables();

        Indexer::Options options;
        options.max_spatial_index_cells(8);
        options.max_spatial_resolution(10 * DEG);
        options.min_spatial_resolution(1 * ARCMIN);
        options.temporal_resolution(10);

        sbs.reindex(options);
        sbs.add_observations(observations);
        sbsdb::add::moving_target(sbs.db(), encke);
        sbsdb::add::observatory(sbs.db(), ztf);
        sbsdb::add::observatory(sbs.db(), rubin);
        sbsdb::add::observatory(sbs.db(), ldt);
    }

    SBSearch<sbsdb::Postgresql> sbs{"postgres:///sbsearch_test"};
    Observations observations{{Observation("test source", "I41", "a", 59252.01, 59252.019, "1:3, 2:3, 2:4, 1:4"),
                               Observation("test source", "I41", "b", 59252.02, 59252.029, "2:3, 3:3, 3:4, 2:4")}};
    sbsearch::MovingTarget encke{"2P"};
    const sbsearch::Observatory ztf{243.14022, 0.836325, +0.546877};
    const sbsearch::Observatory rubin{289.25058, 0.864981, -0.500958};
    const sbsearch::Observatory ldt{248.57749, 0.822887, +0.566916};
};

namespace testing
{
    TEST_F(SBSearchTest, Reindex)
    {
        auto options(sbs.indexer_options());
        options.temporal_resolution(1);
        sbs.reindex(options);

        Observations observations2 = sbsdb::get::observations(sbs.db(), {1, 2});
        EXPECT_NE(observations[0].terms(), observations2[0].terms());
        EXPECT_NE(observations[1].terms(), observations2[1].terms());
    }

    TEST_F(SBSearchTest, FindObservations)
    {
        S2Point point;
        S2Polygon polygon;
        Observations matches;

        // point is observed
        point = S2LatLng::FromDegrees(3.5, 1.5).ToPoint();
        matches = sbs.find_observations(point);
        EXPECT_EQ(matches.size(), 1);

        // and within the time period
        matches = sbs.find_observations(point, {.mjd_start = 59252, .mjd_stop = 59252.021});
        EXPECT_EQ(matches.size(), 1);

        // point is observed, but not within the time period
        matches = sbs.find_observations(point, {.mjd_start = 59252.02});
        EXPECT_EQ(matches.size(), 0);
        matches = sbs.find_observations(point, {.mjd_stop = 59252});
        EXPECT_EQ(matches.size(), 0);

        // point is never observed
        point = S2LatLng::FromDegrees(4.001, 1.5).ToPoint();
        matches = sbs.find_observations(point);
        EXPECT_EQ(matches.size(), 0);

        // invalid time range
        EXPECT_THROW(sbs.find_observations(point, {.mjd_start = 59252.01, .mjd_stop = 59252.00}), std::runtime_error);

        // does not overlap in space
        make_polygon("0:0, 0:1, 1:1", polygon);
        matches = sbs.find_observations(polygon);
        EXPECT_EQ(matches.size(), 0);

        // // but overlaps if padding is given
        matches = sbs.find_observations(polygon, {.padding = 130});
        EXPECT_EQ(matches.size(), 1);

        // does not overlap in space or time
        make_polygon("0:0, 0:1, 1:1", polygon);
        matches = sbs.find_observations(polygon, {.mjd_start = 59252.03, .mjd_stop = 59252.035});
        EXPECT_EQ(matches.size(), 0);

        // overlaps one observation in space
        make_polygon("1:2, 1.5:3.5, 2:2", polygon);
        matches = sbs.find_observations(polygon);
        EXPECT_EQ(matches.size(), 1);

        // overlaps one observation in space, but not time
        make_polygon("1:2, 1.5:3.5, 2:2", polygon);
        matches = sbs.find_observations(polygon, {.mjd_start = 59252.025, .mjd_stop = 59252.035});
        EXPECT_EQ(matches.size(), 0);

        // overlaps one observation in space and time
        make_polygon("1:2, 1.5:3.5, 2:2", polygon);
        matches = sbs.find_observations(polygon, {.mjd_start = 59252.01, .mjd_stop = 59252.022});
        EXPECT_EQ(matches.size(), 1);

        // overlaps two observations in space
        make_polygon("1.5:3, 2.5:3, 2:4", polygon);
        matches = sbs.find_observations(polygon);
        EXPECT_EQ(matches.size(), 2);

        // overlaps two observations in space, but not time
        make_polygon("1.5:3, 2.5:3, 2:4", polygon);
        matches = sbs.find_observations(polygon, {.mjd_start = 59252.05, .mjd_stop = 59252.06});
        EXPECT_EQ(matches.size(), 0);

        // overlaps two observations in space, but only one in time
        make_polygon("1.5:3, 2.5:3, 2:4", polygon);
        matches = sbs.find_observations(polygon, {.mjd_start = 59252.01, .mjd_stop = 59252.022});
        EXPECT_EQ(matches.size(), 1);

        // overlaps two observations in space, and time
        make_polygon("1.5:3, 2.5:3, 2:4", polygon);
        matches = sbs.find_observations(polygon, {.mjd_start = 59252.01, .mjd_stop = 59252.042});
        EXPECT_EQ(matches.size(), 2);

        // invalid time range
        EXPECT_THROW(sbs.find_observations(polygon, {.mjd_start = 59252.01, .mjd_stop = 59252.00}), std::runtime_error);

        // find observations with ephemerides
        // test 1: matches space, but not time
        Ephemeris eph(encke, {{59253.01, 10.01, 0, 3.5, 0, 0, 0, 1, 1, 0},
                              {59253.02, 10.02, 1.5, 3.5, 0, 0, 0, 1, 1, 0},
                              {59253.03, 10.03, 2.5, 3.5, 0, 0, 0, 1, 1, 0},
                              {59253.04, 10.04, 3.5, 3.5, 0, 0, 0, 1, 1, 0}});

        Founds found = sbs.find_observations(eph);
        EXPECT_EQ(found.size(), 0);

        // test 2: matches space and time
        eph = Ephemeris(encke, {{59252.01, 10.01, 0, 3.5, 0, 0, 0, 1, 1, 0},
                                {59252.02, 10.02, 1.5, 3.5, 0, 0, 0, 1, 1, 0},
                                {59252.03, 10.03, 2.5, 3.5, 0, 0, 0, 1, 1, 0},
                                {59252.04, 10.04, 3.5, 3.5, 0, 0, 0, 1, 1, 0}});

        found = sbs.find_observations(eph);
        EXPECT_EQ(found.size(), 2);

        // Add a new data source and limit search by source.
        Observations new_observations({Observation("another test source", "G37", "a", 59252.01, 59252.019, "1:3, 2:3, 2:4, 1:4"),
                                       Observation("another test source", "G37", "b", 59252.02, 59252.029, "2:3, 3:3, 3:4, 2:4")});
        sbs.add_observations(new_observations);
        found = sbs.find_observations(eph);
        EXPECT_EQ(found.size(), 4);

        found = sbs.find_observations(eph, {.source = "test source"});
        EXPECT_EQ(found.size(), 2);

        found = sbs.find_observations(eph, {.source = "another test source"});
        EXPECT_EQ(found.size(), 2);

        // and time
        found = sbs.find_observations(eph, {.mjd_start = 59252.02, .source = "another test source"});
        EXPECT_EQ(found.size(), 1);

        // Search with parallax
        // New observation immediately to the north of previous data
        new_observations = Observations({Observation("another test source", "X05", "c", 59252.01, 59252.019, "1:4, 2:4, 2:5, 1:5")});
        sbs.add_observations(new_observations);
        // this ephemeris is just a few arcsec into the southern FOVs
        eph = Ephemeris(encke, {{59252.01, 10.01, 0.0, 4 - 3.0 / 3600, 0, 0, 0, 1, 1, 0},
                                {59252.02, 10.02, 1.5, 4 - 3.0 / 3600, 0, 0, 0, 1, 1, 0},
                                {59252.03, 10.03, 2.5, 4 - 3.0 / 3600, 0, 0, 0, 1, 1, 0},
                                {59252.04, 10.04, 3.5, 4 - 3.0 / 3600, 0, 0, 0, 1, 1, 0}});

        // nominal geocentric search: expect just the southern FOVs
        found = sbs.find_observations(eph);
        EXPECT_EQ(found.size(), 4);
        auto contains_X05 = [](const Found &f)
        { return f.observation.observatory() == "X05"; };
        EXPECT_EQ(std::count_if(found.begin(), found.end(), contains_X05), 0);

        // add parallax: expect detection in all FOVs
        found = sbs.find_observations(eph, {.parallax = true});
        EXPECT_EQ(found.size(), 5);
        EXPECT_EQ(std::count_if(found.begin(), found.end(), contains_X05), 1);
    }

}
