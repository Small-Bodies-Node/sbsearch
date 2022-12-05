#include "observation.h"
#include "util.h"
#include "sbsearch_testing.h"

#include <cmath>
#include <string>
#include <exception>
#include <s2/s2latlng.h>
#include <s2/s2latlng_rect.h>
#include "s2/s2metrics.h"
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

        TEST(ObservationTests, ObservationIndexTerms)
        {
            // S2LatLngRect rect(S2LatLng::FromDegrees(0, 0), S2LatLng::FromDegrees(1, 1));
            Observation obs(0, 2, "0:0, 0:1, 1:1, 1:0");

            S2RegionTermIndexer::Options options;
            options.set_max_level(S2::kAvgEdge.GetClosestLevel(0.0006)); // 2 deg
            S2RegionTermIndexer indexer(options);

            vector<string> terms = obs.index_terms(indexer);
            // these terms were not independently generated
            for (auto term : sbsearch::split("0555554-0 055555-0 055554-0 05555-0 05554-0 0555-0 0554-0 055-0 $0ffe4-0 0ffe4-0 0fff-0 0ffc-0 0ff-0 $0fffd-0 0fffd-0 0fffc-0 $0ffff-0 0ffff-0 $1001-0 1001-0 1004-0 101-0 $1aaa4-0 1aaa4-0 1aab-0 1aac-0 1ab-0 $1aaa9-0 1aaa9-0 1aaac-0 $1aaab-0 1aaab-0 0555554-1 055555-1 055554-1 05555-1 05554-1 0555-1 0554-1 055-1 $0ffe4-1 0ffe4-1 0fff-1 0ffc-1 0ff-1 $0fffd-1 0fffd-1 0fffc-1 $0ffff-1 0ffff-1 $1001-1 1001-1 1004-1 101-1 $1aaa4-1 1aaa4-1 1aab-1 1aac-1 1ab-1 $1aaa9-1 1aaa9-1 1aaac-1 $1aaab-1 1aaab-1", ' '))
            {
                EXPECT_NE(std::find(terms.begin(), terms.end(), term), terms.end());
            }
        }

        TEST(ObservationTests, ObservationAsPolygonTest)
        {
            Observation obs(0, 1, "-1:-2,2:-2,2:2,-1:2");
            auto polygon = obs.as_polygon();
            auto expected = sbsearch::makePolygon("-1:-2,2:-2,2:2,-1:2");
            EXPECT_TRUE(polygon->Equals(*expected));
        }
    }
}
