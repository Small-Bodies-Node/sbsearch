#include "config.h"

#include <cstdio>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>
#include <gtest/gtest.h>
#include <s2/s2latlng.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>

#include "ephemeris.h"
#include "exceptions.h"
#include "indexer.h"
#include "observation.h"
#include "observatory.h"
#include "sbsearch.h"
#include "util.h"

using sbsearch::Indexer;
using sbsearch::MovingTarget;
using sbsearch::Observation;
using sbsearch::SBSearch;
using std::vector;

class SBSearchTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        Indexer::Options options;
        options.max_spatial_cells(8);
        options.max_spatial_resolution(10 * DEG);
        options.min_spatial_resolution(1 * ARCMIN);
        options.temporal_resolution(10);
        sbs = new SBSearch(SBSearch::sqlite3, ":memory:");
        sbs->reindex(options);
        sbs->add_observations(observations);
        sbs->add_moving_target(encke);
        sbs->add_observatory("I41", ztf);
        sbs->add_observatory("X05", rubin);
        sbs->add_observatory("G37", ldt);
    }

    SBSearch *sbs;
    vector<Observation> observations = {
        Observation("test source", "I41", "a", 59252.01, 59252.019, "1:3, 2:3, 2:4, 1:4"),
        Observation("test source", "I41", "b", 59252.02, 59252.029, "2:3, 3:3, 3:4, 2:4")};
    sbsearch::MovingTarget encke{"2P"};
    const sbsearch::Observatory ztf{243.14022, 0.836325, +0.546877};
    const sbsearch::Observatory rubin{289.25058, 0.864981, -0.500958};
    const sbsearch::Observatory ldt{248.57749, 0.822887, +0.566916};
};

namespace sbsearch
{
    namespace testing
    {
        TEST_F(SBSearchTest, SBSearchStreamInsertOperatorFound)
        {
            Ephemeris eph(encke, {{59252.01, 10.01, 0, 3.5, 0, 0, 0, 1, 1, 0},
                                  {59252.02, 10.02, 1.5, 3.5, 0, 0, 0, 1, 1, 0},
                                  {59252.03, 10.03, 2.5, 3.5, 0, 0, 0, 1, 1, 0},
                                  {59252.04, 10.04, 3.5, 3.5, 0, 0, 0, 1, 1, 0}});
            vector<Found> founds = sbs->find_observations(eph);

            // Should be two found observations
            EXPECT_EQ(founds.size(), 2);

            std::stringstream stream;
            std::string s;
            stream << founds[0];
            EXPECT_EQ(stream.str(), "1  test source  I41  a  59252.01000  59252.01900  777.6  2P  1  59252.01450      0.675000      3.500296   1.000   1.000    0.00");

            stream.str("");
            stream << founds;
            EXPECT_EQ(stream.str(), "observation_id       source  observatory      product_id    mjd_start     mjd_stop  exposure_time  desg  object_id          mjd            ra           dec      rh   delta   phase\n"
                                    "--------------  -----------  -----------  --------------  -----------  -----------  -------------  ----  ---------  -----------  ------------  ------------  ------  ------  ------\n"
                                    "             1  test source          I41               a  59252.01000  59252.01900          777.6    2P          1  59252.01450      0.675000      3.500296   1.000   1.000    0.00\n"
                                    "             2  test source          I41               b  59252.02000  59252.02900          777.6    2P          1  59252.02450      1.950000      3.500132   1.000   1.000    0.00\n");

            stream.str("");
            founds[0].observation.format.show_fov = true;
            stream << founds[0];
            EXPECT_EQ(stream.str(), "1  test source  I41  a  59252.01000  59252.01900  777.6  1:3, 2:3, 2:4, 1:4  2P  1  59252.01450      0.675000      3.500296   1.000   1.000    0.00");

            stream.str("");
            stream << founds;
            EXPECT_EQ(stream.str(), "observation_id       source  observatory      product_id    mjd_start     mjd_stop  exposure_time                 fov  desg  object_id          mjd            ra           dec      rh   delta   phase\n"
                                    "--------------  -----------  -----------  --------------  -----------  -----------  -------------  ------------------  ----  ---------  -----------  ------------  ------------  ------  ------  ------\n"
                                    "             1  test source          I41               a  59252.01000  59252.01900          777.6  1:3, 2:3, 2:4, 1:4    2P          1  59252.01450      0.675000      3.500296   1.000   1.000    0.00\n"
                                    "             2  test source          I41               b  59252.02000  59252.02900          777.6  2:3, 3:3, 3:4, 2:4    2P          1  59252.02450      1.950000      3.500132   1.000   1.000    0.00\n");
        }

        TEST(SBSearchTests, SBSearchReindex)
        {
            char *filename = strdup("/tmp/tmpfileXXXXXX");
            int fd = mkstemp(filename);
            close(fd);

            Indexer::Options options;
            options.max_spatial_cells(8);
            options.max_spatial_resolution(10 * DEG);
            options.min_spatial_resolution(1 * ARCMIN);
            options.temporal_resolution(10);
            SBSearch sbs1(SBSearch::sqlite3, filename);
            sbs1.reindex(options);

            vector<Observation> observations1 = {
                Observation("test source", "I41", "a", 59252.01, 59252.019, "1:3, 2:3, 2:4, 1:4"),
                Observation("test source", "I41", "b", 59252.02, 59252.029, "2:3, 3:3, 3:4, 2:4")};
            sbs1.add_observations(observations1);

            options.temporal_resolution(1);
            sbs1.reindex(options);

            vector<Observation> observations2 = sbs1.get_observations({1, 2});
            EXPECT_NE(observations1[0].terms(), observations2[0].terms());
            EXPECT_NE(observations1[1].terms(), observations2[1].terms());

            std::remove(filename);
        }

        TEST_F(SBSearchTest, SBSearchDateRange)
        {
            vector<Observation> observations{Observation("another test source", "I41", "a", 59253.02, 59253.029, "2:3, 3:3, 3:4, 2:4")};
            sbs->add_observations(observations);

            auto range = sbs->date_range();
            EXPECT_EQ(*range.first, 59252.01);
            EXPECT_EQ(*range.second, 59253.029);

            range = sbs->date_range("test source");
            EXPECT_EQ(*range.first, 59252.01);
            EXPECT_EQ(*range.second, 59252.029);

            range = sbs->date_range("another test source");
            EXPECT_EQ(*range.first, 59253.02);
            EXPECT_EQ(*range.second, 59253.029);

            // nullptr for no observations
            range = sbs->date_range("a third test source");
            EXPECT_EQ(range.first, nullptr);
            EXPECT_EQ(range.second, nullptr);
        }

        TEST_F(SBSearchTest, SBSearchAddGetMovingTargets)
        {
            EXPECT_THROW(sbs->add_moving_target(encke), MovingTargetError);
            MovingTarget halley{"1P"};
            sbs->add_moving_target(halley);

            EXPECT_EQ(sbs->get_moving_target("1P"), halley);
            EXPECT_NE(sbs->get_moving_target("1P"), encke);
            EXPECT_EQ(sbs->get_moving_target("2P"), encke);
            EXPECT_NE(sbs->get_moving_target("2P"), halley);
            EXPECT_EQ(sbs->get_moving_target("2P"), sbs->get_moving_target(1));

            sbs->remove_moving_target(halley);
            EXPECT_THROW(sbs->get_moving_target("1P"), MovingTargetError);
        }

        TEST_F(SBSearchTest, SBSearchAddGetObservatory)
        {
            EXPECT_NO_THROW(sbs->add_observatory("T05", {203.74299, 0.936235, +0.351547}));

            Observatory obs = sbs->get_observatory("I41");
            EXPECT_EQ(obs, ztf);

            Observatories observatories = sbs->get_observatories();
            EXPECT_EQ(observatories["I41"], ztf);

            sbs->remove_observatory("I41");
            EXPECT_THROW(sbs->get_observatory("I41"), ObservatoryError);
        }

        TEST_F(SBSearchTest, SBSearchAddGetEphemerides)
        {
            Ephemeris eph(encke, {{59252.01, 10.01, 0, 3.5, 0, 0, 0, 1, 1, 0},
                                  {59252.02, 10.02, 1.5, 3.5, 0, 0, 0, 1, 1, 0},
                                  {59252.03, 10.03, 2.5, 3.5, 0, 0, 0, 1, 1, 0},
                                  {59252.04, 10.04, 3.5, 3.5, 0, 0, 0, 1, 1, 0}});
            sbs->add_ephemeris(eph);

            Ephemeris test;
            test = sbs->get_ephemeris(encke);
            EXPECT_EQ(test, eph);

            sbs->remove_ephemeris(encke, 59252.025);
            test = sbs->get_ephemeris(encke);
            EXPECT_EQ(test, eph.slice(0, 2));

            test = sbs->get_ephemeris(encke, 0, 59252.015);
            EXPECT_EQ(test, eph[0]);

            sbs->remove_ephemeris(encke, 0, 59252.015);
            test = sbs->get_ephemeris(encke);
            EXPECT_EQ(test, eph[1]);

            sbs->remove_ephemeris(encke);
            test = sbs->get_ephemeris(encke);
            EXPECT_EQ(test.num_vertices(), 0);
        }

        TEST_F(SBSearchTest, SBSearchFindObservations)
        {
            S2Point point;
            S2Polygon polygon;
            vector<Observation> matches;

            // point is observed
            point = S2LatLng::FromDegrees(3.5, 1.5).ToPoint();
            matches = sbs->find_observations(point);
            EXPECT_EQ(matches.size(), 1);

            // and within the time period
            matches = sbs->find_observations(point, {.mjd_start = 59252, .mjd_stop = 59252.021});
            EXPECT_EQ(matches.size(), 1);

            // point is observed, but not within the time period
            matches = sbs->find_observations(point, {.mjd_start = 59252.02});
            EXPECT_EQ(matches.size(), 0);
            matches = sbs->find_observations(point, {.mjd_stop = 59252});
            EXPECT_EQ(matches.size(), 0);

            // point is never observed
            point = S2LatLng::FromDegrees(4.001, 1.5).ToPoint();
            matches = sbs->find_observations(point);
            EXPECT_EQ(matches.size(), 0);

            // invalid time range
            EXPECT_THROW(sbs->find_observations(point, {.mjd_start = 59252.01, .mjd_stop = 59252.00}), std::runtime_error);

            // does not overlap in space
            makePolygon("0:0, 0:1, 1:1", polygon);
            matches = sbs->find_observations(polygon);
            EXPECT_EQ(matches.size(), 0);

            // does not overlap in space or time
            makePolygon("0:0, 0:1, 1:1", polygon);
            matches = sbs->find_observations(polygon, {.mjd_start = 59252.03, .mjd_stop = 59252.035});
            EXPECT_EQ(matches.size(), 0);

            // overlaps one observation in space
            makePolygon("1:2, 1.5:3.5, 2:2", polygon);
            matches = sbs->find_observations(polygon);
            EXPECT_EQ(matches.size(), 1);

            // overlaps one observation in space, but not time
            makePolygon("1:2, 1.5:3.5, 2:2", polygon);
            matches = sbs->find_observations(polygon, {.mjd_start = 59252.025, .mjd_stop = 59252.035});
            EXPECT_EQ(matches.size(), 0);

            // overlaps one observation in space and time
            makePolygon("1:2, 1.5:3.5, 2:2", polygon);
            matches = sbs->find_observations(polygon, {.mjd_start = 59252.01, .mjd_stop = 59252.022});
            EXPECT_EQ(matches.size(), 1);

            // overlaps two observations in space
            makePolygon("1.5:3, 2.5:3, 2:4", polygon);
            matches = sbs->find_observations(polygon);
            EXPECT_EQ(matches.size(), 2);

            // overlaps two observations in space, but not time
            makePolygon("1.5:3, 2.5:3, 2:4", polygon);
            matches = sbs->find_observations(polygon, {.mjd_start = 59252.05, .mjd_stop = 59252.06});
            EXPECT_EQ(matches.size(), 0);

            // overlaps two observations in space, but only one in time
            makePolygon("1.5:3, 2.5:3, 2:4", polygon);
            matches = sbs->find_observations(polygon, {.mjd_start = 59252.01, .mjd_stop = 59252.022});
            EXPECT_EQ(matches.size(), 1);

            // overlaps two observations in space, and time
            makePolygon("1.5:3, 2.5:3, 2:4", polygon);
            matches = sbs->find_observations(polygon, {.mjd_start = 59252.01, .mjd_stop = 59252.042});
            EXPECT_EQ(matches.size(), 2);

            // invalid time range
            EXPECT_THROW(sbs->find_observations(polygon, {.mjd_start = 59252.01, .mjd_stop = 59252.00}), std::runtime_error);

            // find observations with ephemerides
            // test 1: matches space, but not time
            Ephemeris eph(encke, {{59253.01, 10.01, 0, 3.5, 0, 0, 0, 1, 1, 0},
                                  {59253.02, 10.02, 1.5, 3.5, 0, 0, 0, 1, 1, 0},
                                  {59253.03, 10.03, 2.5, 3.5, 0, 0, 0, 1, 1, 0},
                                  {59253.04, 10.04, 3.5, 3.5, 0, 0, 0, 1, 1, 0}});

            vector<Found> found = sbs->find_observations(eph);
            EXPECT_EQ(found.size(), 0);

            // test 2: matches space and time
            eph = Ephemeris(encke, {{59252.01, 10.01, 0, 3.5, 0, 0, 0, 1, 1, 0},
                                    {59252.02, 10.02, 1.5, 3.5, 0, 0, 0, 1, 1, 0},
                                    {59252.03, 10.03, 2.5, 3.5, 0, 0, 0, 1, 1, 0},
                                    {59252.04, 10.04, 3.5, 3.5, 0, 0, 0, 1, 1, 0}});

            found = sbs->find_observations(eph);
            EXPECT_EQ(found.size(), 2);

            // Add a new data source and limit search by source.
            vector<Observation> new_observations{
                Observation("another test source", "G37", "a", 59252.01, 59252.019, "1:3, 2:3, 2:4, 1:4"),
                Observation("another test source", "G37", "b", 59252.02, 59252.029, "2:3, 3:3, 3:4, 2:4")};
            sbs->add_observations(new_observations);
            found = sbs->find_observations(eph);
            EXPECT_EQ(found.size(), 4);

            found = sbs->find_observations(eph, {.source = "test source"});
            EXPECT_EQ(found.size(), 2);

            found = sbs->find_observations(eph, {.source = "another test source"});
            EXPECT_EQ(found.size(), 2);

            // and time
            found = sbs->find_observations(eph, {.mjd_start = 59252.02, .source = "another test source"});
            EXPECT_EQ(found.size(), 1);

            // Search with parallax
            // New observation immediately to the north of previous data
            new_observations = {Observation("another test source", "X05", "c", 59252.01, 59252.019, "1:4, 2:4, 2:5, 1:5")};
            sbs->add_observations(new_observations);
            // this ephemeris is just a few arcsec into the southern FOVs
            eph = Ephemeris(encke, {{59252.01, 10.01, 0.0, 4 - 3.0 / 3600, 0, 0, 0, 1, 1, 0},
                                    {59252.02, 10.02, 1.5, 4 - 3.0 / 3600, 0, 0, 0, 1, 1, 0},
                                    {59252.03, 10.03, 2.5, 4 - 3.0 / 3600, 0, 0, 0, 1, 1, 0},
                                    {59252.04, 10.04, 3.5, 4 - 3.0 / 3600, 0, 0, 0, 1, 1, 0}});
            Observations obs = sbs->find_observations(0, 100000);
            obs[0].format.show_fov = true;

            // nominal geocentric search: expect just the southern FOVs
            found = sbs->find_observations(eph);
            EXPECT_EQ(found.size(), 4);
            auto contains_X05 = [](const Found &f)
            { return f.observation.observatory() == "X05"; };
            EXPECT_EQ(std::count_if(found.begin(), found.end(), contains_X05), 0);

            // add parallax: expect detection in all FOVs
            SBSearch::Options search_options = {.parallax = true, .observatories = sbs->get_observatories()};
            found = sbs->find_observations(eph, search_options);
            EXPECT_EQ(found.size(), 5);
            EXPECT_EQ(std::count_if(found.begin(), found.end(), contains_X05), 1);
        }
    }
}