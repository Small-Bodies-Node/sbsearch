#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <s2/s2polygon.h>

#include "ephemeris.h"
#include "indexer.h"
#include "moving_target.h"
#include "observation.h"
#include "util.h"

using sbsearch::Indexer;
using std::string;
using std::vector;

class IndexerTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        Indexer::Options options;
        options.max_spatial_index_cells(30);
        options.max_spatial_query_cells(8);
        options.max_spatial_resolution(0.17);
        options.min_spatial_resolution(0.01);
        options.temporal_resolution(100);
        indexer = Indexer(options);
    }

    Indexer indexer;
    sbsearch::MovingTarget encke{"2P", 1};
};

namespace sbsearch
{
    namespace testing
    {
        TEST(IndexerTests, OptionsInit)
        {
            Indexer::Options options;

            // check the default values
            EXPECT_EQ(options.max_spatial_index_cells(), 30);
            EXPECT_EQ(options.max_spatial_query_cells(), 8);
            EXPECT_EQ(options.max_spatial_level(), 12);
            EXPECT_EQ(options.min_spatial_level(), 4);
            EXPECT_EQ(options.temporal_resolution(), 1);

            // verify that options are correctly set
            options.max_spatial_index_cells(100);
            EXPECT_EQ(options.max_spatial_index_cells(), 100);

            options.max_spatial_query_cells(10);
            EXPECT_EQ(options.max_spatial_query_cells(), 10);

            options.max_spatial_resolution(10 * ARCMIN);
            EXPECT_NEAR(options.max_spatial_resolution(), 0.002850, 1e-6);
            EXPECT_EQ(options.min_spatial_level(), 9);

            options.min_spatial_resolution(0.1 * ARCMIN);
            EXPECT_NEAR(options.min_spatial_resolution(), 2.227e-05, 1e-8);
            EXPECT_EQ(options.max_spatial_level(), 16);

            options.temporal_resolution(10);
            EXPECT_EQ(options.temporal_resolution(), 10);
        }

        TEST(IndexerTests, OptionsEquality)
        {
            Indexer::Options options;
            Indexer::Options other;

            EXPECT_EQ(options, other);

            other.max_spatial_index_cells(options.max_spatial_index_cells() + 1);
            EXPECT_NE(options, other);

            other.max_spatial_index_cells(options.max_spatial_index_cells());
            EXPECT_EQ(options, other);
            other.max_spatial_query_cells(options.max_spatial_query_cells() + 1);
            EXPECT_NE(options, other);

            other.max_spatial_query_cells(options.max_spatial_query_cells());
            EXPECT_EQ(options, other);
            other.max_spatial_level(options.max_spatial_level() + 1);
            EXPECT_NE(options, other);

            other.max_spatial_level(options.max_spatial_level());
            EXPECT_EQ(options, other);
            other.min_spatial_level(options.min_spatial_level() + 1);
            EXPECT_NE(options, other);

            other.min_spatial_level(options.min_spatial_level());
            EXPECT_EQ(options, other);
            other.temporal_resolution(options.temporal_resolution() + 1);
            EXPECT_NE(options, other);
        }

        TEST(IndexerTests, IndexerInit)
        {
            Indexer indexer;

            // verify default options are set and may be accessed
            EXPECT_EQ(indexer.options().max_spatial_index_cells(), 30);
            EXPECT_EQ(indexer.options().max_spatial_query_cells(), 8);
            EXPECT_EQ(indexer.options().max_spatial_level(), 12);
            EXPECT_EQ(indexer.options().min_spatial_level(), 4);
            EXPECT_EQ(indexer.options().temporal_resolution(), 1);

            Indexer::Options options;
            options.max_spatial_index_cells(100);
            options.max_spatial_query_cells(12);
            options.max_spatial_resolution(0.1);
            options.min_spatial_resolution(0.01);
            options.temporal_resolution(10);
            indexer = Indexer(options);
            EXPECT_EQ(indexer.options().max_spatial_index_cells(), 100);
            EXPECT_EQ(indexer.options().max_spatial_query_cells(), 12);
            EXPECT_EQ(indexer.options().max_spatial_level(), 7);
            EXPECT_EQ(indexer.options().min_spatial_level(), 4);
            EXPECT_EQ(indexer.options().temporal_resolution(), 10);
        }

        TEST_F(IndexerTest, IndexerIndexTermsPoint)
        {
            S2Point point = S2LatLng::FromDegrees(1, 1).ToPoint();

            vector<string> expected = {"1001", "10014", "1004", "101", "104"};

            vector<string> terms = indexer.index_terms(point);
            EXPECT_EQ(std::set<string>(terms.begin(), terms.end()),
                      std::set<string>(expected.begin(), expected.end()));
        }

        TEST_F(IndexerTest, IndexerQueryTermsPoint)
        {
            S2Point point = S2LatLng::FromDegrees(1, 1).ToPoint();

            vector<string> expected = {"$1001", "$10014", "$1004", "$101", "$104", "10014"};

            vector<string> terms = indexer.query_terms(point);
            EXPECT_EQ(std::set<string>(terms.begin(), terms.end()),
                      std::set<string>(expected.begin(), expected.end()));
        }

        TEST_F(IndexerTest, IndexerIndexTermsPolygon)
        {
            S2Polygon polygon;
            make_polygon("1:3, 2:3, 2:4, 1:4", polygon);

            // spatial only
            vector<string> expected = {
                "10194",
                "1019",
                "101c",
                "101",
                "104",
                "1019c",
                "$101b",
                "101b",
                "101c4",
                "101d",
                "101cc",
                "101ec",
                "101f",
            };

            vector<string> terms = indexer.index_terms(polygon);
            EXPECT_EQ(std::set<string>(terms.begin(), terms.end()),
                      std::set<string>(expected.begin(), expected.end()));

            // spatial-temporal
            expected = {
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

            // Here, only expect the first 13 terms
            terms = indexer.index_terms(polygon, 0, 0.01);
            EXPECT_EQ(std::set<string>(terms.begin(), terms.end()),
                      std::set<string>(expected.begin(), expected.begin() + 13));

            // Change the dates and expect all terms
            terms = indexer.index_terms(polygon, 0, 0.02);
            EXPECT_EQ(std::set<string>(terms.begin(), terms.end()),
                      std::set<string>(expected.begin(), expected.end()));
        }

        TEST_F(IndexerTest, IndexerQueryTermsPolygon)
        {
            S2Polygon polygon;
            make_polygon("1:3, 2:3, 2:4, 1:4", polygon);

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

            // Here, only expect the first 13 terms
            vector<string> terms = indexer.query_terms(polygon, 0, 0.01);
            EXPECT_EQ(std::set<string>(terms.begin(), terms.end()),
                      std::set<string>(expected.begin(), expected.begin() + 12));

            // Change the dates and expect all terms
            terms = indexer.query_terms(polygon, 0, 0.02);
            EXPECT_EQ(std::set<string>(terms.begin(), terms.end()),
                      std::set<string>(expected.begin(), expected.end()));
        }

        TEST_F(IndexerTest, IndexerIndexTermsObservation)
        {
            Observation obs("test source", "X05", "product", 0, 0.02, "1:3, 2:3, 2:4, 1:4");
            vector<string> terms = indexer.index_terms(obs);

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

            EXPECT_EQ(std::set<string>(terms.begin(), terms.end()),
                      std::set<string>(expected.begin(), expected.end()));
        }

        TEST_F(IndexerTest, IndexerQueryTermsObservation)
        {
            Observation obs("test source", "X05", "product", 0, 0.02, "1:3, 2:3, 2:4, 1:4");
            vector<string> terms = indexer.query_terms(obs);

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

            EXPECT_EQ(std::set<string>(terms.begin(), terms.end()),
                      std::set<string>(expected.begin(), expected.end()));
        }

        TEST_F(IndexerTest, IndexerIndexTermsEphemeris)
        {
            Ephemeris eph(encke, {{0, 10, 1, 3, 10, 10, 0, 0, 1, 180}, {0.01, 10.01, 2, 4, 10, 10, 0, 2, 1, 90}});
            vector<string> terms = indexer.index_terms(eph);
            std::set<string> expected{
                "101-0",
                "1019-0",
                "10194-0",
                "1019c-0",
                "101b-0",
                "101bc-0",
                "101c-0",
                "101c4-0",
                "101cc-0",
                "101d-0",
                "104-0",
            };
            EXPECT_EQ(std::set<string>(terms.begin(), terms.end()), expected);
        }

        TEST_F(IndexerTest, IndexerQueryTermsEphemeris)
        {
            Ephemeris eph(encke, {{0, 10, 1, 3, 10, 10, 0, 0, 1, 180}, {0.01, 10.01, 2, 4, 10, 10, 0, 2, 1, 90}});
            vector<string> terms = indexer.query_terms(eph);
            std::set<string> expected{
                "$101-0",
                "$1019-0",
                "$101b-0",
                "$101c-0",
                "$101d-0",
                "$104-0",
                "10194-0",
                "1019c-0",
                "101bc-0",
                "101c4-0",
                "101cc-0",
            };
            EXPECT_EQ(std::set<string>(terms.begin(), terms.end()), expected);
        }

    }
}
