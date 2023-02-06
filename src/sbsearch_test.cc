#include "config.h"

#include <vector>
#include <gtest/gtest.h>
#include <s2/s2latlng.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>

#include "ephemeris.h"
#include "indexer.h"
#include "observation.h"
#include "sbsearch.h"
#include "util.h"

using sbsearch::Indexer;
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
        sbs = new SBSearch(SBSearch::sqlite3, ":memory:", options);
        sbs->add_observations(observations);
    }

    SBSearch *sbs;
    vector<Observation> observations = {
        Observation("test source", "a", 59252.01, 59252.019, "1:3, 2:3, 2:4, 1:4"),
        Observation("test source", "b", 59252.02, 59252.029, "2:3, 3:3, 3:4, 2:4")};
};

namespace sbsearch
{
    namespace testing
    {
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
            matches = sbs->find_observations(point, 59252, 59252.018);
            EXPECT_EQ(matches.size(), 1);

            // point is observed, but not within the time period
            matches = sbs->find_observations(point, 59252.02);
            EXPECT_EQ(matches.size(), 0);
            matches = sbs->find_observations(point, -1, 59252);
            EXPECT_EQ(matches.size(), 0);

            // point is never observed
            point = S2LatLng::FromDegrees(4.001, 1.5).ToPoint();
            matches = sbs->find_observations(point);
            EXPECT_EQ(matches.size(), 0);

            // invalid time range
            EXPECT_THROW(sbs->find_observations(point, 59252.01, 59252.00), std::runtime_error);

            // does not overlap in space
            makePolygon("0:0, 0:1, 1:1", polygon);
            matches = sbs->find_observations(polygon);
            EXPECT_EQ(matches.size(), 0);

            // does not overlap in space or time
            makePolygon("0:0, 0:1, 1:1", polygon);
            matches = sbs->find_observations(polygon, 59252.03, 59252.035);
            EXPECT_EQ(matches.size(), 0);

            // overlaps one observation in space
            makePolygon("1:2, 1.5:3.5, 2:2", polygon);
            matches = sbs->find_observations(polygon);
            EXPECT_EQ(matches.size(), 1);

            // overlaps one observation in space, but not time
            makePolygon("1:2, 1.5:3.5, 2:2", polygon);
            matches = sbs->find_observations(polygon, 59252.025, 59252.035);
            EXPECT_EQ(matches.size(), 0);

            // overlaps one observation in space and time
            makePolygon("1:2, 1.5:3.5, 2:2", polygon);
            matches = sbs->find_observations(polygon, 59252.01, 59252.012);
            EXPECT_EQ(matches.size(), 1);

            // overlaps two observations in space
            makePolygon("1.5:3, 2.5:3, 2:4", polygon);
            matches = sbs->find_observations(polygon);
            EXPECT_EQ(matches.size(), 2);

            // overlaps two observations in space, but not time
            makePolygon("1.5:3, 2.5:3, 2:4", polygon);
            matches = sbs->find_observations(polygon, 59252.05, 59252.06);
            EXPECT_EQ(matches.size(), 0);

            // overlaps two observations in space, but only one in time
            makePolygon("1.5:3, 2.5:3, 2:4", polygon);
            matches = sbs->find_observations(polygon, 59252.01, 59252.012);
            EXPECT_EQ(matches.size(), 1);

            // overlaps two observations in space, and time
            makePolygon("1.5:3, 2.5:3, 2:4", polygon);
            matches = sbs->find_observations(polygon, 59252.01, 59252.042);
            EXPECT_EQ(matches.size(), 2);

            // invalid time range
            EXPECT_THROW(sbs->find_observations(polygon, 59252.01, 59252.00), std::runtime_error);

            // find observations with ephemerides
            // test 1: matches space, but not time
            vector<S2Point> vertices{
                S2LatLng::FromDegrees(3.5, 0).ToPoint(), S2LatLng::FromDegrees(3.5, 1.5).ToPoint(),
                S2LatLng::FromDegrees(3.5, 2.5).ToPoint(), S2LatLng::FromDegrees(3.5, 3.5).ToPoint()};
            Ephemeris eph(vertices, {59253.01, 59253.02, 59253.03, 59253.04},
                          {1, 1, 1, 1}, {1, 1, 1, 1}, {0, 0, 0, 0});
            vector<Found> found = sbs->find_observations(eph);
            EXPECT_EQ(found.size(), 0);

            // test 2: matches space and time
            eph = Ephemeris(vertices, {59252.01, 59252.02, 59252.03, 59252.04},
                            {1, 1, 1, 1}, {1, 1, 1, 1}, {0, 0, 0, 0});

            found = sbs->find_observations(eph);
            EXPECT_EQ(found.size(), 2);
        }
    }
}