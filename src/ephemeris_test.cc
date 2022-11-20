#include "ephemeris.h"
#include "util.h"
#include "sbsearch_testing.h"

#include <cmath>
#include <string>
#include <s2/s2point.h>
#include <s2/s2latlng.h>
#include <gtest/gtest.h>

using sbsearch::Ephemeris;
using std::vector;

namespace sbsearch
{
    namespace testing
    {
        TEST(EphemerisTests, EphemerisSegmentTest)
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
        }

        TEST(EphemerisTests, EphemerisTimeTest)
        {
            vector<S2Point> vertices{
                S2LatLng::FromDegrees(0, 1).ToPoint(),
                S2LatLng::FromDegrees(0, 2).ToPoint(),
                S2LatLng::FromDegrees(0, 3).ToPoint()};
            vector<double> times{0, 1, 2};
            EXPECT_NO_THROW(Ephemeris(vertices, times));

            vector<double> badTimes{2, 1, 0};
            EXPECT_ANY_THROW(Ephemeris(vertices, badTimes));
        }

        TEST(EphemerisTests, EphemerisInterpolateTest)
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

        TEST(EphemerisTests, EphemerisInterpolateErrorTest)
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

        TEST(EphemerisTests, EphemerisExtrapolateTest)
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

        TEST(EphemerisTests, EphemerisSubsampleTest)
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

        TEST(EphemerisTests, EphemerisPadTest)
        {
            vector<S2Point> vertices{
                S2LatLng::FromDegrees(0, 0).ToPoint(),
                S2LatLng::FromDegrees(0, 1).ToPoint()};
            vector<double> times{0, 1};
            Ephemeris eph{vertices, times};

            auto polygon = eph.pad(1 * DEG, 2 * DEG);
            auto expected = sbsearch::makePolygon(std::string("-2:-1,-2:0,-2:1,-2:2,2:2,2:1,2:0,2:-1"));
            EXPECT_TRUE(polygon->Equals(*expected));
        }
    }
}