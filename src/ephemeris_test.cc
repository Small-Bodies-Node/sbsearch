#include "config.h"

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

#include "ephemeris.h"
#include "sbsearch_testing.h"
#include "util.h"

using sbsearch::Ephemeris;
using std::vector;

class EphemerisTest : public ::testing::Test
{
protected:
    vector<S2Point> vertices{
        S2LatLng::FromDegrees(0, 1).ToPoint(),
        S2LatLng::FromDegrees(0, 2).ToPoint(),
        S2LatLng::FromDegrees(0, 3).ToPoint()};
    vector<double> mjd{0, 1, 2};
    vector<double> rh{0, 1, 2};
    vector<double> delta{1, 0, 1};
    vector<double> phase{180, 0, 90};
    vector<double> unc_a{1, 5, 10};
    vector<double> unc_b{0.1, 0.5, 1.0};
    vector<double> unc_theta{90, 90, 90};
};

namespace sbsearch
{
    namespace testing
    {
        TEST_F(EphemerisTest, EphemerisInit)
        {
            Ephemeris eph;
            EXPECT_EQ(eph.num_segments(), 0);
            EXPECT_EQ(eph.num_vertices(), 0);

            // single point ephemeris, number of segments should be 0
            eph = Ephemeris(1, {vertices[0]}, {1}, {1}, {1}, {1}, {1}, {1}, {1});
            EXPECT_EQ(eph.num_segments(), 0);
            EXPECT_EQ(eph.num_vertices(), 1);

            eph = Ephemeris(1, vertices, mjd, rh, delta, phase, unc_a, unc_b, unc_theta);
            EXPECT_EQ(vertices, eph.vertices());
            EXPECT_EQ(mjd, eph.mjd());
            EXPECT_EQ(rh, eph.rh());
            EXPECT_EQ(delta, eph.delta());
            EXPECT_EQ(phase, eph.phase());
            EXPECT_EQ(unc_a, eph.unc_a());
            EXPECT_EQ(unc_b, eph.unc_b());
            EXPECT_EQ(unc_theta, eph.unc_theta());

            eph = Ephemeris(1, vertices, mjd, rh, delta, phase);
            vector<double> undefined(3, UNDEF_UNC);
            EXPECT_EQ(eph.unc_a(), undefined);
            EXPECT_EQ(eph.unc_b(), undefined);
            EXPECT_EQ(eph.unc_theta(), undefined);

            // initialize with an invalid mjd array
            EXPECT_THROW(Ephemeris(1, vertices, {3, 2, 1}, rh, delta, phase, unc_a, unc_b, unc_theta), std::runtime_error);

            // initialize with invalid vector sizes
            EXPECT_THROW(Ephemeris(1, {vertices[0], vertices[1]}, mjd, rh, delta, phase, unc_a, unc_b, unc_theta), std::runtime_error);
            EXPECT_THROW(Ephemeris(1, vertices, {0, 1}, rh, delta, phase, unc_a, unc_b, unc_theta), std::runtime_error);
            EXPECT_THROW(Ephemeris(1, vertices, mjd, {0, 1}, delta, phase, unc_a, unc_b, unc_theta), std::runtime_error);
            EXPECT_THROW(Ephemeris(1, vertices, mjd, rh, {0, 1}, phase, unc_a, unc_b, unc_theta), std::runtime_error);
            EXPECT_THROW(Ephemeris(1, vertices, mjd, rh, delta, {0, 1, 2, 3}, unc_a, unc_b, unc_theta), std::runtime_error);
            EXPECT_THROW(Ephemeris(1, vertices, mjd, rh, delta, phase, {0}, unc_b, unc_theta), std::runtime_error);
            EXPECT_THROW(Ephemeris(1, vertices, mjd, rh, delta, phase, unc_a, {}, unc_theta), std::runtime_error);
            EXPECT_THROW(Ephemeris(1, vertices, mjd, rh, delta, phase, unc_a, unc_b, {1, 2, 3, 4, 5, 6}), std::runtime_error);
        }

        TEST_F(EphemerisTest, EphemerisVertex)
        {
            Ephemeris eph = Ephemeris(1, vertices, mjd, rh, delta, phase, unc_a, unc_b, unc_theta);
            for (int i = 0; i < vertices.size(); i++)
                EXPECT_EQ(eph.vertex(i), eph.vertex(i - 3));

            EXPECT_THROW(eph.vertex(3), std::runtime_error);
            EXPECT_THROW(eph.vertex(-4), std::runtime_error);
        }

        TEST_F(EphemerisTest, EphemerisPropertyGetter)
        {
            Ephemeris eph = Ephemeris(1, vertices, mjd, rh, delta, phase, unc_a, unc_b, unc_theta);
            EXPECT_THROW(eph.mjd(3), std::runtime_error);
        }

        TEST_F(EphemerisTest, EphemerisIsEqual)
        {
            Ephemeris eph(1, vertices, mjd, rh, delta, phase, unc_a, unc_b, unc_theta);
            Ephemeris same(1, vertices, mjd, rh, delta, phase, unc_a, unc_b, unc_theta);
            EXPECT_TRUE(eph.is_equal(same));
            EXPECT_TRUE(same.is_equal(eph));

            Ephemeris not_same = eph.segment(0);
            EXPECT_FALSE(eph.is_equal(not_same));
            EXPECT_FALSE(not_same.is_equal(eph));

            vector<S2Point> different_vertices{
                S2LatLng::FromDegrees(0, 1).ToPoint(),
                S2LatLng::FromDegrees(0, 1.5).ToPoint(),
                S2LatLng::FromDegrees(0, 3).ToPoint()};
            not_same = Ephemeris(1, different_vertices, mjd, rh, delta, phase, unc_a, unc_b, unc_theta);
            EXPECT_FALSE(eph.is_equal(not_same));
            EXPECT_FALSE(not_same.is_equal(eph));

            not_same = Ephemeris(1, vertices, {0, 0.5, 2}, rh, delta, phase, unc_a, unc_b, unc_theta);
            EXPECT_FALSE(eph.is_equal(not_same));
            EXPECT_FALSE(not_same.is_equal(eph));

            not_same = Ephemeris(1, vertices, mjd, rh, {0, 0.5, 2}, phase, unc_a, unc_b, unc_theta);
            EXPECT_FALSE(eph.is_equal(not_same));
            EXPECT_FALSE(not_same.is_equal(eph));

            not_same = Ephemeris(1, vertices, mjd, rh, delta, {0, 0.5, 2}, unc_a, unc_b, unc_theta);
            EXPECT_FALSE(eph.is_equal(not_same));
            EXPECT_FALSE(not_same.is_equal(eph));

            not_same = Ephemeris(1, vertices, mjd, rh, delta, phase, {0, 0.5, 2}, unc_b, unc_theta);
            EXPECT_FALSE(eph.is_equal(not_same));
            EXPECT_FALSE(not_same.is_equal(eph));

            not_same = Ephemeris(1, vertices, mjd, rh, delta, phase, unc_a, {0, 0.5, 2}, unc_theta);
            EXPECT_FALSE(eph.is_equal(not_same));
            EXPECT_FALSE(not_same.is_equal(eph));

            not_same = Ephemeris(1, vertices, mjd, rh, delta, phase, unc_a, unc_b, {0, 0.5, 2});
            EXPECT_FALSE(eph.is_equal(not_same));
            EXPECT_FALSE(not_same.is_equal(eph));
        }

        TEST_F(EphemerisTest, EphemerisBracketOperator)
        {
            Ephemeris eph(1, vertices, mjd, rh, delta, phase, unc_a, unc_b, unc_theta);

            for (int i = 0; i < eph.num_vertices(); i++)
            {
                Ephemeris a = eph[i];
                EXPECT_EQ(a.vertex(0), eph.vertex(i));
                EXPECT_EQ(a.mjd(0), eph.mjd(i));
                EXPECT_EQ(a.rh(0), eph.rh(i));
                EXPECT_EQ(a.delta(0), eph.delta(i));
                EXPECT_EQ(a.phase(0), eph.phase(i));
                EXPECT_EQ(a.unc_a(0), eph.unc_a(i));
                EXPECT_EQ(a.unc_b(0), eph.unc_b(i));
                EXPECT_EQ(a.unc_theta(0), eph.unc_theta(i));

                a = eph[i - eph.num_vertices()];
                EXPECT_EQ(a.vertex(0), eph.vertex(i));
                EXPECT_EQ(a.mjd(0), eph.mjd(i));
                EXPECT_EQ(a.rh(0), eph.rh(i));
                EXPECT_EQ(a.delta(0), eph.delta(i));
                EXPECT_EQ(a.phase(0), eph.phase(i));
                EXPECT_EQ(a.unc_a(0), eph.unc_a(i));
                EXPECT_EQ(a.unc_b(0), eph.unc_b(i));
                EXPECT_EQ(a.unc_theta(0), eph.unc_theta(i));
            }

            EXPECT_THROW(eph[eph.num_vertices()], std::runtime_error);
            EXPECT_THROW(eph[-eph.num_vertices() - 1], std::runtime_error);
        }

        TEST_F(EphemerisTest, EphemerisAppend)
        {
            Ephemeris eph(1, vertices, mjd, rh, delta, phase, unc_a, unc_b, unc_theta);
            Ephemeris a(1, {S2LatLng::FromDegrees(1, 4).ToPoint()},
                        {3}, {3}, {3}, {45});
            EXPECT_EQ(eph.object_id(), a.object_id());
            eph.append(a);
            EXPECT_TRUE(eph[3].is_equal(a));

            // append to an empty ephemeris
            Ephemeris b;
            EXPECT_THROW(b.append(a), std::runtime_error);
            b.object_id(a.object_id()); // need the object ids to match
            b.append(a);
            EXPECT_TRUE(a.is_equal(b));

            // don't append if the mjd is out of order
            EXPECT_THROW(a.append(eph), std::runtime_error);
        }

        TEST_F(EphemerisTest, EphemerisSegment)
        {
            Ephemeris eph(1, vertices, mjd, rh, delta, phase, unc_a, unc_b, unc_theta);

            // get single segment
            Ephemeris segment = eph.segment(1);
            Ephemeris expected(1,
                               {vertices[1], vertices[2]},
                               {mjd[1], mjd[2]},
                               {rh[1], rh[2]},
                               {delta[1], delta[2]},
                               {phase[1], phase[2]},
                               {unc_a[1], unc_a[2]},
                               {unc_b[1], unc_b[2]},
                               {unc_theta[1], unc_theta[2]});
            EXPECT_TRUE(segment.is_equal(expected));

            // get all segments
            vector<Ephemeris> segments = eph.segments();
            for (int i = 0; i < eph.num_segments(); i++)
            {
                Ephemeris s = eph.segment(i);
                EXPECT_TRUE(segments[i].is_equal(s));

                s = eph.segment(i - eph.num_segments());
                EXPECT_TRUE(segments[i].is_equal(s));
            }

            // fail on invalid index
            EXPECT_THROW(eph.segment(eph.num_segments()), std::runtime_error);
            EXPECT_THROW(eph.segment(-eph.num_segments() - 1), std::runtime_error);
        }

        TEST_F(EphemerisTest, EphemerisAsPolyline)
        {
            Ephemeris eph(1, vertices, mjd, rh, delta, phase, unc_a, unc_b, unc_theta);
            S2Polyline polyline(vertices);
            EXPECT_TRUE(eph.as_polyline().Equals(polyline));
        }

        TEST_F(EphemerisTest, EphemerisInterpolate)
        {
            Ephemeris eph(1, vertices, mjd, rh, delta, phase, unc_a, unc_b, unc_theta);

            Ephemeris a = eph.interpolate(0.5);
            S2LatLng coord = S2LatLng(a.vertex(0));
            EXPECT_NEAR(coord.lat().degrees(), 0, 1 * ARCSEC);
            EXPECT_NEAR(coord.lng().degrees(), 1.5, 1 * ARCSEC);
            EXPECT_EQ(a.mjd(0), 0.5);
            EXPECT_EQ(a.rh(0), 0.5);
            EXPECT_EQ(a.delta(0), 0.5);
            EXPECT_EQ(a.phase(0), 90);
            EXPECT_EQ(a.unc_a(0), 3);
            EXPECT_NEAR(a.unc_b(0), 0.3, 1e-8);
            EXPECT_EQ(a.unc_theta(0), 90);

            a = eph.interpolate(1.5);
            coord = S2LatLng(a.vertex(0));
            EXPECT_NEAR(coord.lat().degrees(), 0, 1 * ARCSEC);
            EXPECT_NEAR(coord.lng().degrees(), 2.5, 1 * ARCSEC);

            EXPECT_THROW(eph.interpolate(-1), std::runtime_error);
            EXPECT_THROW(eph.interpolate(3), std::runtime_error);
        }

        TEST_F(EphemerisTest, EphemerisExtrapolate)
        {
            using Extrapolate = Ephemeris::Extrapolate;
            Ephemeris eph(1, vertices, mjd, rh, delta, phase, unc_a, unc_b, unc_theta);

            Ephemeris a = eph.extrapolate(5 * DEG, Extrapolate::BACKWARDS);
            S2LatLng coord = S2LatLng(a.vertex(0));
            EXPECT_NEAR(coord.lat().degrees(), 0, 1 * ARCSEC);
            EXPECT_NEAR(coord.lng().degrees(), -4, 1 * ARCSEC);

            a = eph.extrapolate(5 * DEG, Extrapolate::FORWARDS);
            coord = S2LatLng(a.vertex(0));
            EXPECT_NEAR(coord.lat().degrees(), 0, 1 * ARCSEC);
            EXPECT_NEAR(coord.lng().degrees(), 8, 1 * ARCSEC);
        }

        TEST_F(EphemerisTest, EphemerisSubsample)
        {
            Ephemeris eph(1, vertices, mjd, rh, delta, phase, unc_a, unc_b, unc_theta);

            Ephemeris subsample = eph.subsample(0.5, 0.75);
            EXPECT_EQ(subsample.num_segments(), 1);
            EXPECT_TRUE(subsample[0].is_equal(eph.interpolate(0.5)));
            EXPECT_TRUE(subsample[1].is_equal(eph.interpolate(0.75)));

            subsample = eph.subsample(0.5, 1.5);
            EXPECT_EQ(subsample.num_segments(), 2);
            EXPECT_TRUE(subsample[0].is_equal(eph.interpolate(0.5)));
            EXPECT_TRUE(subsample[1].is_equal(eph[1]));
            EXPECT_TRUE(subsample[2].is_equal(eph.interpolate(1.5)));

            subsample = eph.subsample(1, 2);
            EXPECT_EQ(subsample.num_segments(), 1);
            EXPECT_TRUE(subsample[0].is_equal(eph[1]));
            EXPECT_TRUE(subsample[1].is_equal(eph[2]));
        }

        S2Polygon generate_expected_polygon(const S2Point &start, const S2Point &end, const double a, const double b, const double theta)
        {
            // a, b, theta in radians

            // get the ellipses for the first and last points
            vector<S2LatLng> e0 = ellipse(16, S2LatLng(start), a, b, theta);
            vector<S2LatLng> e1 = ellipse(16, S2LatLng(end), a, b, theta);

            // ephmeris vector is along RA, and our padded region is elongated along Dec
            vector<S2LatLng> coords;
            // start with half of the first ellipse
            coords.insert(coords.end(), e0.begin() + 8, e0.end());
            coords.push_back(e0[0]);
            // append half the last ellipse
            coords.insert(coords.end(), e1.begin(), e1.begin() + 9);

            vector<S2Point> points(coords.size());
            std::transform(coords.begin(), coords.end(), points.begin(), [](S2LatLng c)
                           { return c.ToPoint(); });

            S2Polygon expected;
            makePolygon(points, expected);
            return expected;
        }

        TEST_F(EphemerisTest, EphemerisPad)
        {
            Ephemeris eph(1, vertices, mjd, rh, delta, phase, unc_a, unc_b, unc_theta);

            auto polygon = eph.pad(3600 * 2, 3600, 0);

            // which is equivalent to:
            auto polygon2 = eph.pad(3600, 2 * 3600);
            EXPECT_TRUE(polygon.BoundaryNear(polygon2, S1Angle::Radians(1 * ARCSEC)));

            S2Polygon expected = generate_expected_polygon(vertices[0], vertices[2], 2 * DEG, 1 * DEG, 0);
            EXPECT_TRUE(polygon.BoundaryNear(expected, S1Angle::Radians(1 * ARCSEC)));

            EXPECT_THROW(eph.pad({1}, {1, 2, 3}), std::runtime_error);
            EXPECT_THROW(eph.pad({1, 2, 3}, {1, 2}, {0, 0, 0}), std::runtime_error);
            EXPECT_THROW(eph.pad(3600 * 90, 3600 * 90), std::runtime_error);
        }

        TEST_F(EphemerisTest, EphemerisAsPolygon)
        {
            Ephemeris eph(1, vertices, mjd, rh, delta, phase, {10, 10, 10}, {10, 10, 10}, unc_theta);

            eph.mutable_options()->padding = 3600;
            S2Polygon polygon = eph.as_polygon();

            S2Polygon expected_polygon = generate_expected_polygon(vertices[0], vertices[2], 1 * DEG, 1 * DEG, 0);
            EXPECT_TRUE(polygon.BoundaryNear(expected_polygon, S1Angle::Radians(1 * ARCSEC)));

            eph.mutable_options()->padding = 0;
            eph.mutable_options()->use_uncertainty = true;
            polygon = eph.as_polygon();
            expected_polygon = generate_expected_polygon(vertices[0], vertices[2], 10 * ARCSEC, 10 * ARCSEC, 0);
            EXPECT_TRUE(polygon.BoundaryNear(expected_polygon, S1Angle::Radians(1 * ARCSEC)));
        }
    }
}