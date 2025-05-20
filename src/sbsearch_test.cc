#include "config.h"

#include <string>
#include <vector>
#include <gtest/gtest.h>
#include <s2/s2latlng.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>

#include "constants.h"
#include "ephemeris.h"
#include "exceptions.h"
#include "indexer.h"
#include "observation.h"
#include "sbsdb/sbsdb.h"
#include "sbsearch.h"
#include "util/polygon.h"
#include "util/string.h"

using namespace sbsearch;
using sbsearch::util::make_polygon;
using sbsearch::util::make_vertices;
using std::endl;
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

    SBSearch<sbsdb::Postgresql> sbs{"postgres:///sbsearch_test", {.log_level = LogLevel::DEBUG}};
    Observations observations{{Observation("test source", "I41", "a", 59252.01, 59252.019, "1:3, 2:3, 2:4, 1:4"),
                               Observation("test source", "I41", "b", 59252.02, 59252.029, "2:3, 3:3, 3:4, 2:4")}};
    sbsearch::MovingTarget encke{"2P"};
    const sbsearch::Observatory ztf{243.14022, 0.836325, +0.546877, "I41"};
    const sbsearch::Observatory rubin{289.25058, 0.864981, -0.500958, "X05"};
    const sbsearch::Observatory ldt{248.57749, 0.822887, +0.566916, "G37"};
};

namespace testing
{
    TEST_F(SBSearchTest, Reindex)
    {
        auto options(sbs.indexer_options());
        options.max_spatial_resolution(1 * DEG);
        options.min_spatial_resolution(10 * ARCMIN);
        sbs.reindex(options);

        for (int i = 0; i < 2; i++)
        {
            auto obs = sbsdb::get::observations(sbs.db(), {i + 1})[0];
            // terms should be different
            EXPECT_NE(observations[i].terms(), obs.terms());
            // center should be the same
            EXPECT_EQ(observations[i].center(), obs.center());
        }
    }

    TEST_F(SBSearchTest, AddObservations)
    {
        // verify that the indices were added

        // Checked the center positions with python.  There is no guarantee that
        // the terms will be the same in different versions of s2
        /*
        >>> import pywraps2 as s2
        >>> from astropy.coordinates import SkyCoord
        >>> indexer = s2.S2RegionTermIndexer()
        >>> indexer.set_min_level(30)
        >>> indexer.set_max_level(30)
        >>> c = SkyCoord([1, 2, 2, 1], [3, 3, 4, 4], unit="deg")
        >>> c = SkyCoord(c.cartesian.x.mean(), c.cartesian.y.mean(), c.cartesian.z.mean(), representation_type="cartesian").spherical
        >>> c.lat
        <Latitude 3.50013294 deg>
        >>> c.lon
        <Longitude 1.5 deg>
        >>> point = s2.S2LatLng.FromDegrees(c.lat.deg, c.lon.deg).ToPoint()
        >>> indexer.GetIndexTerms(point, "")
        ('101be98a94063031',)
        >>> c = SkyCoord([2, 3, 3, 2], [3, 3, 4, 4], unit="deg")
        >>> c = SkyCoord(c.cartesian.x.mean(), c.cartesian.y.mean(), c.cartesian.z.mean(), representation_type="cartesian").spherical>>> c.lat
        <Latitude 3.50013294 deg>
        >>> c.lon
        <Longitude 2.5 deg>
        >>> point = s2.S2LatLng.FromDegrees(c.lat.deg, c.lon.deg).ToPoint()
        >>> indexer.GetIndexTerms(point, "")
        ('1010a6d587f7a239',)
        */

        Observations observations = sbsdb::get::observations(sbs.db(), {1, 2});
        EXPECT_EQ(observations[0].terms(), vector<string>({"$10195", "10195", "10194", "1019", "101c", "101", "104", "$10197", "10197", "$10199", "10199", "1019c", "$101b", "101b", "$101c1", "101c1", "101c4", "101d", "$101c7", "101c7", "$101c9", "101c9", "101cc", "$101eb", "101eb", "101ec", "101f"}));
        EXPECT_EQ(observations[0].center(), "101be98a94063031");

        EXPECT_EQ(observations[1].terms(), vector<string>({"$10104", "10104", "1011", "1014", "101", "104", "$1010c", "1010c", "$10173", "10173", "10174", "1017", "$10175", "10175", "$1019c", "1019c", "1019", "101c", "$101a4", "101a4", "101b", "$101a9", "101a9", "101ac", "$101af4", "101af4", "101af"}));
        EXPECT_EQ(observations[1].center(), "1010a6d587f7a239");
    }

    TEST_F(SBSearchTest, FindObservationsByPoint)
    {
        S2Point point;
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
    }

    TEST_F(SBSearchTest, FindObservationsByPolygon)
    {
        S2Polygon polygon;
        Observations matches;

        // does not overlap in space
        make_polygon(make_vertices("0:0, 0:1, 1:1"), polygon);
        matches = sbs.find_observations(polygon);
        EXPECT_EQ(matches.size(), 0);

        // but overlaps if padding is given
        matches = sbs.find_observations(polygon, {.padding = 130});
        EXPECT_EQ(matches.size(), 1);

        // does not overlap in space or time
        make_polygon(make_vertices("0:0, 0:1, 1:1"), polygon);
        matches = sbs.find_observations(polygon, {.mjd_start = 59252.03, .mjd_stop = 59252.035});
        EXPECT_EQ(matches.size(), 0);

        // overlaps one observation in space
        make_polygon(make_vertices("1:2, 1.5:3.5, 2:2"), polygon);
        matches = sbs.find_observations(polygon);
        EXPECT_EQ(matches.size(), 1);

        // overlaps one observation in space, but not time
        make_polygon(make_vertices("1:2, 1.5:3.5, 2:2"), polygon);
        matches = sbs.find_observations(polygon, {.mjd_start = 59252.025, .mjd_stop = 59252.035});
        EXPECT_EQ(matches.size(), 0);

        // overlaps one observation in space and time
        make_polygon(make_vertices("1:2, 1.5:3.5, 2:2"), polygon);
        matches = sbs.find_observations(polygon, {.mjd_start = 59252.01, .mjd_stop = 59252.022});
        EXPECT_EQ(matches.size(), 1);

        // overlaps two observations in space
        make_polygon(make_vertices("1.5:3, 2.5:3, 2:4"), polygon);
        matches = sbs.find_observations(polygon);
        EXPECT_EQ(matches.size(), 2);

        // overlaps two observations in space, but not time
        make_polygon(make_vertices("1.5:3, 2.5:3, 2:4"), polygon);
        matches = sbs.find_observations(polygon, {.mjd_start = 59252.05, .mjd_stop = 59252.06});
        EXPECT_EQ(matches.size(), 0);

        // overlaps two observations in space, but only one in time
        make_polygon(make_vertices("1.5:3, 2.5:3, 2:4"), polygon);
        matches = sbs.find_observations(polygon, {.mjd_start = 59252.01, .mjd_stop = 59252.022});
        EXPECT_EQ(matches.size(), 1);

        // overlaps two observations in space, and time
        make_polygon(make_vertices("1.5:3, 2.5:3, 2:4"), polygon);
        matches = sbs.find_observations(polygon, {.mjd_start = 59252.01, .mjd_stop = 59252.042});
        EXPECT_EQ(matches.size(), 2);

        // invalid time range
        EXPECT_THROW(sbs.find_observations(polygon, {.mjd_start = 59252.01, .mjd_stop = 59252.00}), std::runtime_error);
    }

    TEST_F(SBSearchTest, FindObservationsByEphemeris)
    {
        // find observations with ephemerides
        Logger::debug() << "\n    Test 1: matches space, but not time" << endl;
        Ephemeris eph(encke, {{59253.01, 10.01, 0, 3.5, 0, 0, 0, 1, 1, 0},
                              {59253.02, 10.02, 1.5, 3.5, 0, 0, 0, 1, 1, 0},
                              {59253.03, 10.03, 2.5, 3.5, 0, 0, 0, 1, 1, 0},
                              {59253.04, 10.04, 3.5, 3.5, 0, 0, 0, 1, 1, 0}});

        Founds found = sbs.find_observations(eph);
        EXPECT_EQ(found.size(), 0);

        Logger::debug() << "\n    Test 2: matches space and time" << endl;
        eph = Ephemeris(encke, {{59252.01, 10.01, 0, 3.5, 0, 0, 0, 1, 1, 0},
                                {59252.02, 10.02, 1.5, 3.5, 0, 0, 0, 1, 1, 0},
                                {59252.03, 10.03, 2.5, 3.5, 0, 0, 0, 1, 1, 0},
                                {59252.04, 10.04, 3.5, 3.5, 0, 0, 0, 1, 1, 0}});

        found = sbs.find_observations(eph);
        EXPECT_EQ(found.size(), 2);

        Logger::debug() << "\n    Add a new data source and limit search by source" << endl;
        Observations new_observations({Observation("another test source", "G37", "c", 59252.01, 59252.019, "1:3, 2:3, 2:4, 1:4"),
                                       Observation("another test source", "G37", "d", 59252.02, 59252.029, "2:3, 3:3, 3:4, 2:4")});
        sbs.add_observations(new_observations);
        found = sbs.find_observations(eph);
        EXPECT_EQ(found.size(), 4);

        found = sbs.find_observations(eph, {.source = "test source"});
        EXPECT_EQ(found.size(), 2);

        found = sbs.find_observations(eph, {.source = "another test source"});
        EXPECT_EQ(found.size(), 2);

        Logger::debug() << "\n    Add a new data source and limit search by source and time" << endl;
        found = sbs.find_observations(eph, {.mjd_start = 59252.02, .source = "another test source"});
        EXPECT_EQ(found.size(), 1);

        Logger::debug() << "\n    Search with parallax" << endl;
        // New observation immediately to the north of previous data
        new_observations = Observations({Observation("another test source", "X05", "e", 59252.01, 59252.019, "1:4, 2:4, 2:5, 1:5")});
        sbs.add_observations(new_observations);
        // this ephemeris is just a few arcsec into the southern FOVs
        eph = Ephemeris(encke, {{59252.01, 10.01, 0.0, 4 - 3.0 / 3600, 0, 0, 0, 1, 1, 0},
                                {59252.02, 10.02, 1.5, 4 - 3.0 / 3600, 0, 0, 0, 1, 1, 0},
                                {59252.03, 10.03, 2.5, 4 - 3.0 / 3600, 0, 0, 0, 1, 1, 0},
                                {59252.04, 10.04, 3.5, 4 - 3.0 / 3600, 0, 0, 0, 1, 1, 0}});

        Logger::debug() << "\n    Nominal geocentric search: expect just the southern FOVs" << endl;
        found = sbs.find_observations(eph);
        EXPECT_EQ(found.size(), 4);
        auto contains_X05 = [](const Found &f)
        { return f.observation.observatory() == "X05"; };
        EXPECT_EQ(std::count_if(found.begin(), found.end(), contains_X05), 0);

        Logger::debug() << "\n    Add parallax: expect detection in all FOVs" << endl;
        found = sbs.find_observations(eph, {.parallax = true});
        EXPECT_EQ(found.size(), 5);
        EXPECT_EQ(std::count_if(found.begin(), found.end(), contains_X05), 1);
    }
}
