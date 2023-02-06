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
            Observation obs("test source", "product", 0, 0.1, "0:0, 0:1, 1:1");
            EXPECT_TRUE(obs.is_valid());
            EXPECT_EQ(obs.observation_id(), UNDEFINED_OBSID);
            EXPECT_EQ(obs.mjd_start(), 0);
            EXPECT_EQ(obs.mjd_stop(), 0.1);

            obs = Observation("test source", "product", 1, 1.1, "0:0, 0:1, 1:1", "", 1);
            EXPECT_EQ(obs.observation_id(), 1);

            // bad FOV: only two vertices
            EXPECT_THROW(Observation("test source", "product", 0, 0.1, "0:0, 0:1"), std::runtime_error);
            // bad FOV: not parsable as coordinates
            EXPECT_THROW(Observation("test source", "product", 0, 0.1, "asdf"), std::runtime_error);
            // stop is before start
            EXPECT_THROW(Observation("test source", "product", 0.1, 0, "0:0, 0:1, 1:1"), std::runtime_error);
        }

        TEST(ObservationTests, ObservationInitTestS2LatLng)
        {
            vector<S2LatLng> vertices{
                S2LatLng::FromDegrees(0, 0),
                S2LatLng::FromDegrees(0, 1),
                S2LatLng::FromDegrees(1, 1)};
            Observation obs("test source", "product", 0, 0.1, vertices);
            EXPECT_TRUE(obs.is_valid());
            EXPECT_EQ(obs.fov(), "0.000000:0.000000, 1.000000:0.000000, 1.000000:1.000000");
        }

        TEST(ObservationTests, ObservationProperties)
        {
            Observation a("test source 1", "a", 1, 1.1, "0:0, 0:1, 1:1", "asdf fdsa", 1);
            EXPECT_EQ(a.mjd_start(), 1);
            EXPECT_EQ(a.mjd_stop(), 1.1);
            EXPECT_EQ(a.fov(), "0:0, 0:1, 1:1");
            EXPECT_EQ(a.terms(), "asdf fdsa");

            Observation b("test source 2", "b", 2, 2.1, "2:0, 2:1, 3:1", "jkl; ;lkj", 1);
            EXPECT_FALSE(a == b);

            a.source("test source 2");
            EXPECT_FALSE(a == b);

            a.product_id("b");
            EXPECT_FALSE(a == b);

            a.mjd_start(2);
            EXPECT_FALSE(a == b);

            a.mjd_stop(2.1);
            EXPECT_FALSE(a == b);

            a.fov("2:0, 2:1, 3:1");
            EXPECT_TRUE(a == b); // terms are not compared

            // cannot update observation_id
            EXPECT_THROW(a.observation_id(2), std::runtime_error);

            // but can update it if it was undefined
            Observation c("test source", "c", 1, 1.1, "0:0, 0:1, 1:1");
            EXPECT_EQ(c.observation_id(), UNDEFINED_OBSID);
            c.observation_id(2);
            EXPECT_EQ(c.observation_id(), 2);

            // update terms from a vector
            a.terms(vector<string>{"a", "b", "c"});
            EXPECT_EQ(a.terms(), "a b c");
        }

        TEST(ObservationTests, ObservationTerms)
        {
            vector<S2LatLng> vertices{
                S2LatLng::FromDegrees(0, 0),
                S2LatLng::FromDegrees(0, 1),
                S2LatLng::FromDegrees(1, 1)};
            Observation obs("test source", "obs", 0, 0.1, vertices, "asdf fsda");
            EXPECT_EQ(obs.terms(), "asdf fsda");
        }

        TEST(ObservationTests, ObservationIsSameFov)
        {
            vector<S2LatLng> vertices{
                S2LatLng::FromDegrees(0, 0),
                S2LatLng::FromDegrees(1, 0),
                S2LatLng::FromDegrees(1, 1)};
            Observation a("test source", "a", 0, 0.1, vertices);
            Observation b("test source", "b", 1, 1.1, "0:0, 0:1, 1:1");
            EXPECT_FALSE(a.is_equal(b));
            EXPECT_TRUE(a.is_same_fov(b));
        }

        TEST(ObservationTests, ObservationIsEqual)
        {
            vector<S2LatLng> vertices{
                S2LatLng::FromDegrees(0, 0),
                S2LatLng::FromDegrees(1, 0),
                S2LatLng::FromDegrees(1, 1)};
            Observation obs("test source", "a", 0, 0.1, vertices, "", 1);

            // same
            Observation other("test source", "a", 0, 0.1, "0:0, 0:1, 1:1", "", 1);
            EXPECT_TRUE(obs.is_equal(other));

            // different source
            other = Observation("another test source", "a", 0, 0.1, "0:0, 0:1, 1:1", "", 1);
            EXPECT_FALSE(obs.is_equal(other));

            // different product_id
            other = Observation("test source", "b", 0, 0.1, "0:0, 0:1, 1:1", "", 1);
            EXPECT_FALSE(obs.is_equal(other));

            // different mjd_start
            other = Observation("test source", "a", 0.05, 0.1, "0:0, 0:1, 1:1", "", 1);
            EXPECT_FALSE(obs.is_equal(other));

            // different mjd_stop
            other = Observation("test source", "a", 0, 0.15, "0:0, 0:1, 1:1", "", 1);
            EXPECT_FALSE(obs.is_equal(other));

            // different observation_id
            other = Observation("test source", "a", 0, 0.1, "0:0, 0:1, 1:1", "", 2);
            EXPECT_FALSE(obs.is_equal(other));

            // different terms
            other = Observation("test source", "a", 0, 0.1, "0:0, 0:1, 1:1", "asdf", 1);
            EXPECT_TRUE(obs.is_equal(other));

            // different FOV
            other = Observation("test source", "a", 0, 0.1, "0.05:0, 0:1, 1:1", "asdf", 1);
            EXPECT_FALSE(obs.is_equal(other));
        }

        TEST(ObservationTests, ObservationEqualityOperator)
        {
            vector<S2LatLng> vertices{
                S2LatLng::FromDegrees(0, 0),
                S2LatLng::FromDegrees(1, 0),
                S2LatLng::FromDegrees(1, 1)};
            Observation obs("test source", "a", 0, 0.1, vertices, "", 1);

            // same
            Observation other("test source", "a", 0, 0.1, "0:0, 0:1, 1:1", "", 1);
            EXPECT_TRUE(obs == other);
        }

        TEST(ObservationTests, ObservationObservationId)
        {
            Observation obs("test source", "product", 0, 0.1, "0:0, 0:1, 1:1");
            obs.observation_id(1);
            EXPECT_EQ(obs.observation_id(), 1);

            obs = Observation("test source", "product", 1, 1.1, "0:0, 0:1, 1:1", "", 1);
            EXPECT_THROW(obs.observation_id(2), std::runtime_error);
        }

        TEST(ObservationTests, ObservationAsPolygon)
        {
            Observation obs("test source", "product", 0, 1, "-1:-2,2:-2,2:2,-1:2");
            auto polygon = obs.as_polygon();
            S2Polygon expected;
            makePolygon("-1:-2,2:-2,2:2,-1:2", expected);
            EXPECT_TRUE(polygon.Equals(expected));
        }
    }
}
