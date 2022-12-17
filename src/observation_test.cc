#include "observation.h"
#include "util.h"
#include "sbsearch_testing.h"

#include <cmath>
#include <string>
#include <exception>
#include <s2/s2latlng.h>
#include <s2/s2latlng_rect.h>
#include <s2/s2metrics.h>
#include <s2/s2point.h>
#include <s2/s2region_term_indexer.h>
#include <gtest/gtest.h>

using sbsearch::Observation;
using std::string;
using std::vector;

namespace sbsearch
{
    namespace testing
    {
        TEST(ObservationTests, ObservationInitTestFov)
        {
            Observation obs(0, 0.1, "0:0, 0:1, 1:1");
            EXPECT_TRUE(obs.is_valid());
            EXPECT_EQ(obs.observation_id(), UNDEFINED_OBSID);
            EXPECT_EQ(obs.mjd_start(), 0);
            EXPECT_EQ(obs.mjd_stop(), 0.1);

            obs = Observation(1, 1.1, "0:0, 0:1, 1:1", "", 1);
            EXPECT_EQ(obs.observation_id(), 1);

            EXPECT_THROW(Observation obs(0, 0.1, "0:0, 0:1"), std::runtime_error);
            EXPECT_THROW(Observation obs(0, 0.1, "asdf"), std::runtime_error);
        }

        TEST(ObservationTests, ObservationInitTestS2LatLng)
        {
            vector<S2LatLng> vertices{
                S2LatLng::FromDegrees(0, 0),
                S2LatLng::FromDegrees(0, 1),
                S2LatLng::FromDegrees(1, 1)};
            Observation obs(0, 0.1, vertices);
            EXPECT_TRUE(obs.is_valid());
        }

        TEST(ObservationTests, ObservationFov)
        {
            vector<S2LatLng> vertices{
                S2LatLng::FromDegrees(0, 0),
                S2LatLng::FromDegrees(0, 1),
                S2LatLng::FromDegrees(1, 1)};
            Observation obs(0, 0.1, vertices);
            EXPECT_EQ(obs.fov(), "0.000000:0.000000, 1.000000:0.000000, 1.000000:1.000000");
        }

        TEST(ObservationTests, ObservationTerms)
        {
            vector<S2LatLng> vertices{
                S2LatLng::FromDegrees(0, 0),
                S2LatLng::FromDegrees(0, 1),
                S2LatLng::FromDegrees(1, 1)};
            Observation obs(0, 0.1, vertices, "asdf, fsda");
            EXPECT_EQ(obs.terms(), "asdf, fsda"); // not yet generated
        }

        TEST(ObservationTests, ObservationIsSameFov)
        {
            vector<S2LatLng> vertices{
                S2LatLng::FromDegrees(0, 0),
                S2LatLng::FromDegrees(1, 0),
                S2LatLng::FromDegrees(1, 1)};
            Observation a(0, 0.1, vertices);
            Observation b(1, 1.1, "0:0, 0:1, 1:1");
            EXPECT_TRUE(a.is_same_fov(b));
        }

        TEST(ObservationTests, ObservationIsEqual)
        {
            vector<S2LatLng> vertices{
                S2LatLng::FromDegrees(0, 0),
                S2LatLng::FromDegrees(1, 0),
                S2LatLng::FromDegrees(1, 1)};
            Observation obs(0, 0.1, vertices, "", 1);

            // same
            Observation other(0, 0.1, "0:0, 0:1, 1:1", "", 1);
            EXPECT_TRUE(obs.is_equal(other));

            // different mjd_start
            other = Observation(0.05, 0.1, "0:0, 0:1, 1:1", "", 1);
            EXPECT_FALSE(obs.is_equal(other));

            // different mjd_stop
            other = Observation(0, 0.15, "0:0, 0:1, 1:1", "", 1);
            EXPECT_FALSE(obs.is_equal(other));

            // different observation_id
            other = Observation(0, 0.1, "0:0, 0:1, 1:1", "", 2);
            EXPECT_FALSE(obs.is_equal(other));

            // different terms
            other = Observation(0, 0.1, "0:0, 0:1, 1:1", "asdf", 1);
            EXPECT_TRUE(obs.is_equal(other));

            // different FOV
            other = Observation(0, 0.1, "0.05:0, 0:1, 1:1", "asdf", 1);
            EXPECT_FALSE(obs.is_equal(other));
        }

        TEST(ObservationTests, ObservationObservationId)
        {
            Observation obs(0, 0.1, "0:0, 0:1, 1:1");
            obs.observation_id(1);
            EXPECT_EQ(obs.observation_id(), 1);

            obs = Observation(1, 1.1, "0:0, 0:1, 1:1", "", 1);
            EXPECT_THROW(obs.observation_id(2), std::runtime_error);
        }

        TEST(ObservationTests, ObservationIndexTerms)
        {
            vector<string> expected = {
                "10194-0",
                "1019-0",
                "101c-0",
                "101-0",
                "104-0",
                "1019c-0",
                "$101b-0",
                "101b-0",
                "101c4-0",
                "101d-0",
                "101cc-0",
                "101ec-0",
                "101f-0",
                "10194-1",
                "1019-1",
                "101c-1",
                "101-1",
                "104-1",
                "1019c-1",
                "$101b-1",
                "101b-1",
                "101c4-1",
                "101d-1",
                "101cc-1",
                "101ec-1",
                "101f-1",
            };

            S2RegionTermIndexer::Options options;
            options.set_min_level(S2::kAvgEdge.GetClosestLevel(0.17));
            options.set_max_level(S2::kAvgEdge.GetClosestLevel(0.01));
            options.set_max_cells(8);
            S2RegionTermIndexer indexer(options);

            // Here, only expect the first 13 terms
            Observation obs(0, 0.1, "1:3, 2:3, 2:4, 1:4");
            vector<string> terms = obs.index_terms(indexer);
            EXPECT_EQ(std::set<string>(terms.begin(), terms.end()),
                      std::set<string>(expected.begin(), expected.begin() + 13));

            // Now, expect all terms
            obs = Observation(0, 2, "1:3, 2:3, 2:4, 1:4");
            terms = obs.index_terms(indexer);
            EXPECT_EQ(std::set<string>(terms.begin(), terms.end()),
                      std::set<string>(expected.begin(), expected.end()));
        }

        TEST(ObservationTests, ObservationQueryTerms)
        {
            vector<string> expected = {
                "$101-0",
                "$1019-0",
                "$101c-0",
                "$101d-0",
                "$101f-0",
                "$104-0",
                "10194-0",
                "1019c-0",
                "101b-0",
                "101c4-0",
                "101cc-0",
                "101ec-0",
                "$101-1",
                "$1019-1",
                "$101c-1",
                "$101d-1",
                "$101f-1",
                "$104-1",
                "10194-1",
                "1019c-1",
                "101b-1",
                "101c4-1",
                "101cc-1",
                "101ec-1",
            };

            S2RegionTermIndexer::Options options;
            options.set_min_level(S2::kAvgEdge.GetClosestLevel(0.17));
            options.set_max_level(S2::kAvgEdge.GetClosestLevel(0.01));
            options.set_max_cells(8);
            S2RegionTermIndexer indexer(options);

            // Here, only expect the first 12 terms
            Observation obs(0, 0.1, "1:3, 2:3, 2:4, 1:4");
            vector<string> terms = obs.query_terms(indexer);
            EXPECT_EQ(std::set<string>(terms.begin(), terms.end()),
                      std::set<string>(expected.begin(), expected.begin() + 12));

            // Now, expect all terms
            obs = Observation(0, 2, "1:3, 2:3, 2:4, 1:4");
            terms = obs.query_terms(indexer);
            EXPECT_EQ(std::set<string>(terms.begin(), terms.end()),
                      std::set<string>(expected.begin(), expected.end()));
        }

        TEST(ObservationTests, ObservationAsPolygon)
        {
            Observation obs(0, 1, "-1:-2,2:-2,2:2,-1:2");
            auto polygon = obs.as_polygon();
            auto expected = sbsearch::makePolygon("-1:-2,2:-2,2:2,-1:2");
            EXPECT_TRUE(polygon->Equals(*expected));
        }
    }
}
