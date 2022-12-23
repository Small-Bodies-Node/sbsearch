#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <s2/s2polygon.h>

#include "indexer.h"
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
        options.max_spatial_cells(8);
        options.max_spatial_resolution(0.17);
        options.min_spatial_resolution(0.01);
        indexer = Indexer(options);
    }

    Indexer indexer;
};

namespace sbsearch
{
    namespace testing
    {
        TEST(IndexerTests, OptionsInit)
        {
            Indexer::Options options;

            // check the default values
            EXPECT_EQ(options.max_spatial_cells(), 8);
            EXPECT_EQ(options.max_spatial_level(), 12);
            EXPECT_EQ(options.min_spatial_level(), 4);
            EXPECT_EQ(options.temporal_resolution(), 100);

            // verify that options are correctly set
            options.max_spatial_cells(12);
            EXPECT_EQ(options.max_spatial_cells(), 12);

            options.max_spatial_resolution(10 * ARCMIN);
            EXPECT_NEAR(options.max_spatial_resolution(), 0.002850, 1e-6);
            EXPECT_EQ(options.min_spatial_level(), 9);

            options.min_spatial_resolution(0.1 * ARCMIN);
            EXPECT_NEAR(options.min_spatial_resolution(), 2.227e-05, 1e-8);
            EXPECT_EQ(options.max_spatial_level(), 16);

            options.temporal_resolution(10);
            EXPECT_EQ(options.temporal_resolution(), 10);
        }

        TEST(IndexerTests, IndexerInit)
        {
            Indexer indexer;

            // verify default options are set and may be accessed
            EXPECT_EQ(indexer.options().max_spatial_cells(), 8);
            EXPECT_EQ(indexer.options().max_spatial_level(), 12);
            EXPECT_EQ(indexer.options().min_spatial_level(), 4);
            EXPECT_EQ(indexer.options().temporal_resolution(), 100);

            Indexer::Options options;
            options.max_spatial_cells(8);
            options.max_spatial_resolution(0.1);
            options.min_spatial_resolution(0.01);
            options.temporal_resolution(10);
            indexer = Indexer(options);
            EXPECT_EQ(indexer.options().max_spatial_cells(), 8);
            EXPECT_EQ(indexer.options().max_spatial_level(), 7);
            EXPECT_EQ(indexer.options().min_spatial_level(), 4);
            EXPECT_EQ(indexer.options().temporal_resolution(), 10);
        }

        // {
        //     // test 2 time resolutions, 11 points per time resolution, width of 3 points
        //     // expect: first 8 steps have one term, next 3 have two terms, repeat
        //     for (int step = 0; step < 2 * 11; step++)
        //     {
        //         const double start = 59800.0 + step / 11.0 / TIME_TERMS_PER_DAY;
        //         const double stop = 59800.0 + (step + 3) / 11.0 / TIME_TERMS_PER_DAY;
        //         const vector<string> terms = mjd_to_time_terms(start, stop);

        //         EXPECT_EQ(terms[0], std::to_string(59800 * TIME_TERMS_PER_DAY + (unsigned int)floor(step / 11.0)));
        //         if (step % 11 < 9)
        //         {
        //             EXPECT_EQ(terms.size(), 1);
        //         }
        //         else
        //         {
        //             EXPECT_EQ(terms.size(), 2);
        //             EXPECT_EQ(terms[1], std::to_string(59800 * TIME_TERMS_PER_DAY + (unsigned int)floor((step + 3) / 11.0)));
        //         }
        //     }
        // }

        TEST_F(IndexerTest, IndexerIndexTerms)
        {
            S2Polygon polygon;
            makePolygon("1:3, 2:3, 2:4, 1:4", polygon);

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

            // Here, only expect the first 13 terms
            vector<string> terms = indexer.index_terms(polygon, 0, 0.01);
            EXPECT_EQ(std::set<string>(terms.begin(), terms.end()),
                      std::set<string>(expected.begin(), expected.begin() + 13));

            // Change the dates and expect all terms
            terms = indexer.index_terms(polygon, 0, 0.02);
            EXPECT_EQ(std::set<string>(terms.begin(), terms.end()),
                      std::set<string>(expected.begin(), expected.end()));
        }

        TEST_F(IndexerTest, IndexerQueryTerms)
        {
            S2Polygon polygon;
            makePolygon("1:3, 2:3, 2:4, 1:4", polygon);

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
    }
}
