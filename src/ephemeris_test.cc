#include "ephemeris.h"
#include "sbsearch_testing.h"
#include "util.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <string>
#include <s2/s2latlng.h>
#include <s2/s2metrics.h>
#include <s2/s2point.h>
#include <s2/s2polyline.h>
#include <s2/s2region_term_indexer.h>
#include <gtest/gtest.h>

using sbsearch::Ephemeris;
using std::vector;

namespace sbsearch
{
    namespace testing
    {
        TEST(EphemerisTests, EphemerisIsEqual)
        {
            vector<S2Point> vertices{
                S2LatLng::FromDegrees(0, 1).ToPoint(),
                S2LatLng::FromDegrees(0, 2).ToPoint(),
                S2LatLng::FromDegrees(0, 3).ToPoint()};
            vector<double> times{0, 1, 2};
            Ephemeris eph{vertices, times};
            Ephemeris same{vertices, times};
            EXPECT_TRUE(eph.is_equal(same));
            EXPECT_TRUE(same.is_equal(eph));

            Ephemeris not_same = eph.segment(0);
            EXPECT_FALSE(eph.is_equal(not_same));
            EXPECT_FALSE(not_same.is_equal(eph));

            vector<S2Point> different_vertices{
                S2LatLng::FromDegrees(0, 1).ToPoint(),
                S2LatLng::FromDegrees(0, 1.5).ToPoint(),
                S2LatLng::FromDegrees(0, 3).ToPoint()};
            not_same = Ephemeris(different_vertices, times);
            EXPECT_FALSE(eph.is_equal(not_same));
            EXPECT_FALSE(not_same.is_equal(eph));

            vector<double> different_times{0, 0.5, 2};
            not_same = Ephemeris(vertices, different_times);
            EXPECT_FALSE(eph.is_equal(not_same));
            EXPECT_FALSE(not_same.is_equal(eph));
        }

        TEST(EphemerisTests, EphemerisVertex)
        {
            vector<S2Point> vertices{
                S2LatLng::FromDegrees(0, 1).ToPoint(),
                S2LatLng::FromDegrees(0, 2).ToPoint(),
                S2LatLng::FromDegrees(0, 3).ToPoint()};
            vector<double> times{0, 1, 2};
            Ephemeris eph{vertices, times};
            vector<S2Point> v = eph.vertices();
            EXPECT_EQ(v[0], vertices[0]);
            EXPECT_EQ(v[1], vertices[1]);
            EXPECT_EQ(v[2], vertices[2]);
        }

        TEST(EphemerisTests, EphemerisSegment)
        {
            vector<S2Point> vertices{
                S2LatLng::FromDegrees(0, 1).ToPoint(),
                S2LatLng::FromDegrees(0, 2).ToPoint(),
                S2LatLng::FromDegrees(0, 3).ToPoint()};
            vector<double> times{0, 1, 2};
            Ephemeris eph{vertices, times};
            Ephemeris segment = eph.segment(1);

            vector<S2Point> vertices2{
                S2LatLng::FromDegrees(0, 2).ToPoint(),
                S2LatLng::FromDegrees(0, 3).ToPoint()};
            vector<double> times2{1, 2};
            Ephemeris expected{vertices2, times2};

            EXPECT_EQ(segment.vertex(0), expected.vertex(0));
            EXPECT_EQ(segment.vertex(1), expected.vertex(1));
            EXPECT_EQ(segment.time(0), expected.time(0));
            EXPECT_EQ(segment.time(1), expected.time(1));

            vector<Ephemeris> segments = eph.segments();
            for (int i = 0; i < 2; i++)
            {
                Ephemeris s = eph.segment(i);
                EXPECT_TRUE(segments[i].is_equal(s));
            }
        }

        TEST(EphemerisTests, EphemerisTime)
        {
            vector<S2Point> vertices{
                S2LatLng::FromDegrees(0, 1).ToPoint(),
                S2LatLng::FromDegrees(0, 2).ToPoint(),
                S2LatLng::FromDegrees(0, 3).ToPoint()};
            vector<double> times{0, 1, 2};
            Ephemeris eph(vertices, times);

            vector<double> t = eph.times();
            for (int i = 0; i < 3; i++)
                EXPECT_EQ(t[i], times[i]);

            vector<double> badTimes{2, 1, 0};
            EXPECT_ANY_THROW(Ephemeris(vertices, badTimes));
        }

        TEST(EphemerisTests, EphemerisAsPolyline)
        {
            vector<S2Point> vertices{
                S2LatLng::FromDegrees(0, 1).ToPoint(),
                S2LatLng::FromDegrees(0, 2).ToPoint(),
                S2LatLng::FromDegrees(0, 3).ToPoint()};
            vector<double> times{0, 1, 2};
            Ephemeris eph(vertices, times);

            S2Polyline polyline(vertices);
            EXPECT_TRUE(eph.as_polyline().Equals(polyline));
        }

        TEST(EphemerisTests, EphemerisInterpolate)
        {
            vector<S2Point> vertices{
                S2LatLng::FromDegrees(0, 0).ToPoint(),
                S2LatLng::FromDegrees(0, 20).ToPoint(),
                S2LatLng::FromDegrees(10, 20).ToPoint()};
            vector<double> times{0, 1, 2};
            Ephemeris eph{vertices, times};

            S2LatLng coord = S2LatLng(eph.interpolate(0.5));
            EXPECT_NEAR(coord.lat().degrees(), 0, 1 * ARCSEC);
            EXPECT_NEAR(coord.lng().degrees(), 10, 1 * ARCSEC);

            coord = S2LatLng(eph.interpolate(1.5));
            EXPECT_NEAR(coord.lat().degrees(), 5, 1 * ARCSEC);
            EXPECT_NEAR(coord.lng().degrees(), 20, 1 * ARCSEC);
        }

        TEST(EphemerisTests, EphemerisInterpolateError)
        {
            vector<S2Point> vertices{
                S2LatLng::FromDegrees(0, 0).ToPoint(),
                S2LatLng::FromDegrees(0, 20).ToPoint(),
                S2LatLng::FromDegrees(10, 20).ToPoint()};
            vector<double> times{0, 1, 2};
            Ephemeris eph{vertices, times};

            EXPECT_ANY_THROW(eph.interpolate(-1));
            EXPECT_ANY_THROW(eph.interpolate(3));
        }

        TEST(EphemerisTests, EphemerisExtrapolate)
        {
            using Extrapolate = Ephemeris::Extrapolate;

            vector<S2Point> vertices{
                S2LatLng::FromDegrees(0, 0).ToPoint(),
                S2LatLng::FromDegrees(0, 20).ToPoint(),
                S2LatLng::FromDegrees(10, 20).ToPoint()};
            vector<double> times{0, 1, 2};
            Ephemeris eph{vertices, times};

            S2LatLng coord = S2LatLng(eph.extrapolate(5 * DEG, Extrapolate::BACKWARDS));
            EXPECT_NEAR(coord.lat().degrees(), 0, 1 * ARCSEC);
            EXPECT_NEAR(coord.lng().degrees(), -5, 1 * ARCSEC);

            coord = S2LatLng(eph.extrapolate(5 * DEG, Extrapolate::FORWARDS));
            EXPECT_NEAR(coord.lat().degrees(), 15, 1 * ARCSEC);
            EXPECT_NEAR(coord.lng().degrees(), 20, 1 * ARCSEC);
        }

        TEST(EphemerisTests, EphemerisSubsample)
        {
            vector<S2Point> vertices{
                S2LatLng::FromDegrees(0, 0).ToPoint(),
                S2LatLng::FromDegrees(0, 20).ToPoint(),
                S2LatLng::FromDegrees(10, 20).ToPoint(),
                S2LatLng::FromDegrees(20, 20).ToPoint()};
            vector<double> times{0, 1, 2, 3};
            Ephemeris eph{vertices, times};

            Ephemeris subsample = eph.subsample(0.5, 1.5);
            EXPECT_EQ(subsample.num_segments(), 1);
            EXPECT_EQ(subsample.vertex(0), eph.interpolate(0.5));
            EXPECT_EQ(subsample.vertex(1), eph.interpolate(1.5));

            subsample = eph.subsample(0.5, 2.5);
            EXPECT_EQ(subsample.num_segments(), 3);
            EXPECT_EQ(subsample.vertex(0), eph.interpolate(0.5));
            EXPECT_EQ(subsample.vertex(1), eph.vertex(1));
            EXPECT_EQ(subsample.vertex(2), eph.vertex(2));
            EXPECT_EQ(subsample.vertex(3), eph.interpolate(2.5));

            subsample = eph.subsample(1, 2);
            EXPECT_EQ(subsample.num_segments(), 1);
            EXPECT_EQ(subsample.vertex(0), eph.vertex(1));
            EXPECT_EQ(subsample.vertex(1), eph.vertex(2));

            subsample = eph.subsample(0.25, 0.75);
            EXPECT_EQ(subsample.num_segments(), 1);
            EXPECT_EQ(subsample.vertex(0), eph.interpolate(0.25));
            EXPECT_EQ(subsample.vertex(1), eph.interpolate(0.75));
        }

        TEST(EphemerisTests, EphemerisPad)
        {
            vector<S2Point> vertices{
                S2LatLng::FromDegrees(0, 0).ToPoint(),
                S2LatLng::FromDegrees(0, 1).ToPoint()};
            vector<double> times{0, 1};
            Ephemeris eph{vertices, times};

            auto polygon = eph.pad(1 * DEG, 2 * DEG);
            auto expected = sbsearch::makePolygon(std::string("-1:-2,0:-2,1:-2,2:-2,2:2,1:2,0:2,-1:2"));
            EXPECT_TRUE(polygon->Equals(*expected));
        }

        TEST(EphemerisTests, EphemerisQueryTerms)
        {
            vector<S2Point> vertices{
                S2LatLng::FromDegrees(3, 1).Normalized().ToPoint(),
                S2LatLng::FromDegrees(4, 2).Normalized().ToPoint()};
            vector<double> times{0, 1};
            Ephemeris eph{vertices, times};

            S2RegionTermIndexer::Options options;
            options.set_min_level(S2::kAvgEdge.GetClosestLevel(0.17));
            options.set_max_level(S2::kAvgEdge.GetClosestLevel(0.01));
            options.set_max_cells(8);
            S2RegionTermIndexer indexer(options);

            vector<string> terms = eph.query_terms(indexer);
            // expected values generated with the python code
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

        TEST(EphemerisTests, EphemerisQueryTermsPadded)
        {
            S2RegionTermIndexer::Options options1;
            options1.set_min_level(S2::kAvgEdge.GetClosestLevel(0.17));
            options1.set_max_level(S2::kAvgEdge.GetClosestLevel(0.01));
            options1.set_max_cells(8);
            S2RegionTermIndexer indexer2(options1);

            vector<S2Point> vertices{
                S2LatLng::FromDegrees(0, 0).ToPoint(),
                S2LatLng::FromDegrees(0, 0.01).ToPoint()};
            vector<double> times{0, 1};
            Ephemeris eph{vertices, times};
            vector<string> terms = eph.query_terms(indexer2, 0.01 * DEG, 0.01 * DEG);

            // query terms should match this region
            auto polygon = makePolygon("-0.01:0.01, 0.02:0.01, 0.02:-0.01, -0.01:-0.01");
            vector<string> expected = indexer2.GetQueryTerms(*polygon, "");
            std::transform(
                expected.begin(), expected.end(), expected.begin(),
                [](string s)
                { return s + "-0"; });

            EXPECT_EQ(std::set<string>(terms.begin(), terms.end()),
                      std::set<string>(expected.begin(), expected.end()));
        }

        TEST(EphemerisTests, EphemerisIndexTerms)
        {
            vector<S2Point> vertices{
                S2LatLng::FromDegrees(3, 1).Normalized().ToPoint(),
                S2LatLng::FromDegrees(4, 2).Normalized().ToPoint()};
            vector<double> times{0, 1};
            Ephemeris eph{vertices, times};

            S2RegionTermIndexer::Options options;
            options.set_min_level(S2::kAvgEdge.GetClosestLevel(0.17));
            options.set_max_level(S2::kAvgEdge.GetClosestLevel(0.01));
            options.set_max_cells(8);
            S2RegionTermIndexer indexer(options);

            vector<string> terms = eph.index_terms(indexer);
            // expected values generated with index_terms, and is not a true independent test
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

        TEST(EphemerisTests, EphemerisIndexTermsPadded)
        {
            S2RegionTermIndexer::Options options;
            options.set_min_level(S2::kAvgEdge.GetClosestLevel(0.17));
            options.set_max_level(S2::kAvgEdge.GetClosestLevel(0.01));
            options.set_max_cells(8);
            S2RegionTermIndexer indexer(options);

            vector<S2Point> vertices{
                S2LatLng::FromDegrees(0, 0).ToPoint(),
                S2LatLng::FromDegrees(0, 0.01).ToPoint()};
            vector<double> times{0, 1};
            Ephemeris eph{vertices, times};
            vector<string> terms = eph.index_terms(indexer, 0.01 * DEG, 0.01 * DEG);

            // query terms should match this region
            auto polygon = makePolygon("-0.01:0.01, 0.02:0.01, 0.02:-0.01, -0.01:-0.01");
            vector<string> expected = indexer.GetIndexTerms(*polygon, "");
            std::transform(
                expected.begin(), expected.end(), expected.begin(),
                [](string s)
                { return s + "-0"; });

            EXPECT_EQ(std::set<string>(terms.begin(), terms.end()),
                      std::set<string>(expected.begin(), expected.end()));
        }

    }
}