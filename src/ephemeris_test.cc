#include "config.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <string>
#include <sstream>
#include <s2/s2latlng.h>
#include <s2/s2metrics.h>
#include <s2/s2point.h>
#include <s2/s2polyline.h>
#include <s2/s2region_term_indexer.h>
#include <gtest/gtest.h>

#include "ephemeris.h"
#include "moving_target.h"
#include "sbsearch_testing.h"
#include "util.h"

using sbsearch::Ephemeris;
using std::vector;

class EphemerisTest : public ::testing::Test
{
protected:
    Ephemeris::Data data{
        // mjd, tmtp, ra, dec, unc a, b, theta, rh, delta, phase, selong, true, sangle, vangle, vmag
        {0, 10, 1, 0, 1, 0.1, 90, 0, 1, 180, 0, 0, 0, 10, -1},
        {1, 11, 2, 0, 5, 0.5, 90, 1, 0, 0, 180, 30, 0, 20, 5},
        {2, 12, 3, 0, 10, 1.0, 90, 2, 1, 90, 80, 90, 0, 30, 10}};
    sbsearch::MovingTarget encke{"2P", 1};
};

namespace sbsearch
{
    namespace testing
    {
        TEST(EphemerisTests, EphemerisDatumEquality)
        {
            Ephemeris::Datum a{0, 10, 1, 0, 1, 0.1, 90, 0, 1, 180, 0, 0, 0, 10, -1};
            Ephemeris::Datum b{0, 10, 1, 0, 1, 0.1, 90, 0, 1, 180, 0, 0, 0, 10, -1};
            Ephemeris::Datum c{1, 11, 2, 0, 5, 0.5, 90, 1, 0, 0, 180, 30, 0, 20, 5};
            EXPECT_EQ(a, b);
            EXPECT_NE(a, c);
            EXPECT_NE(b, c);

            Ephemeris::Datum d{0, 10, 1, 0, 1, 0.1, 90, 0, 1, 180, 0, 0, 0, 10, -1};
            EXPECT_EQ(a, d);
            d.mjd += 1;
            EXPECT_NE(a, d);

            d.mjd = a.mjd;
            EXPECT_EQ(a, d);
            d.tmtp += 1;
            EXPECT_NE(a, d);

            d.tmtp = a.tmtp;
            EXPECT_EQ(a, d);
            d.ra += 1;
            EXPECT_NE(a, d);

            d.ra = a.ra;
            EXPECT_EQ(a, d);
            d.dec += 1;
            EXPECT_NE(a, d);

            d.dec = a.dec;
            EXPECT_EQ(a, d);
            d.unc_a += 1;
            EXPECT_NE(a, d);

            d.unc_a = a.unc_a;
            EXPECT_EQ(a, d);
            d.unc_b += 1;
            EXPECT_NE(a, d);

            d.unc_b = a.unc_b;
            EXPECT_EQ(a, d);
            d.unc_theta += 1;
            EXPECT_NE(a, d);

            d.unc_theta = a.unc_theta;
            EXPECT_EQ(a, d);
            d.rh += 1;
            EXPECT_NE(a, d);

            d.rh = a.rh;
            EXPECT_EQ(a, d);
            d.delta += 1;
            EXPECT_NE(a, d);

            d.delta = a.delta;
            EXPECT_EQ(a, d);
            d.phase += 1;
            EXPECT_NE(a, d);

            d.phase = a.phase;
            EXPECT_EQ(a, d);
            d.selong += 1;
            EXPECT_NE(a, d);

            d.selong = a.selong;
            EXPECT_EQ(a, d);
            d.true_anomaly += 1;
            EXPECT_NE(a, d);

            d.true_anomaly = a.true_anomaly;
            EXPECT_EQ(a, d);
            d.sangle += 1;
            EXPECT_NE(a, d);

            d.sangle = a.sangle;
            EXPECT_EQ(a, d);
            d.vangle += 1;
            EXPECT_NE(a, d);

            d.vangle = a.vangle;
            EXPECT_EQ(a, d);
            d.vmag += 1;
            EXPECT_NE(a, d);
        }

        TEST(EphemerisTests, EphemerisDatumRaDecFromPoint)
        {
            Ephemeris::Datum d;
            d.radec(S2LatLng::FromDegrees(0, 1).ToPoint());
            EXPECT_EQ(d.ra, 1);
            EXPECT_EQ(d.dec, 0);
        }

        TEST_F(EphemerisTest, EphemerisDatumAsPoint)
        {
            EXPECT_EQ(data[0].as_s2point(), S2LatLng::FromDegrees(0, 1).ToPoint());
        }

        TEST_F(EphemerisTest, EphemerisDatumAsJSON)
        {
            json::object obj = data[0].as_json();
            EXPECT_EQ(obj["mjd"], 0.);
            EXPECT_EQ(obj["tmtp"], 10.);
            EXPECT_EQ(obj["ra"], 1.);
            EXPECT_EQ(obj["dec"], 0.);
            EXPECT_EQ(obj["unc_a"], 1.);
            EXPECT_EQ(obj["unc_b"], 0.1);
            EXPECT_EQ(obj["unc_theta"], 90.);
            EXPECT_EQ(obj["rh"], 0.);
            EXPECT_EQ(obj["delta"], 1.);
            EXPECT_EQ(obj["phase"], 180.);
            EXPECT_EQ(obj["selong"], 0.);
            EXPECT_EQ(obj["true_anomaly"], 0.);
            EXPECT_EQ(obj["sangle"], 0.);
            EXPECT_EQ(obj["vangle"], 10.);
            EXPECT_EQ(obj["vmag"], -1.);
        }

        TEST_F(EphemerisTest, EphemerisInit)
        {
            Ephemeris eph;
            EXPECT_EQ(eph.num_segments(), 0);
            EXPECT_EQ(eph.num_vertices(), 0);

            // single point ephemeris, number of segments should be 0
            eph = Ephemeris(encke, {data[0]});
            EXPECT_EQ(eph.num_segments(), 0);
            EXPECT_EQ(eph.num_vertices(), 1);

            eph = Ephemeris(encke, data);
            EXPECT_EQ(vector<double>({0, 1, 2}), eph.mjd());
            EXPECT_EQ(vector<double>({10, 11, 12}), eph.tmtp());
            EXPECT_EQ(vector<double>({1, 2, 3}), eph.ra());
            EXPECT_EQ(vector<double>({0, 0, 0}), eph.dec());
            EXPECT_EQ(vector<double>({1, 5, 10}), eph.unc_a());
            EXPECT_EQ(vector<double>({0.1, 0.5, 1.0}), eph.unc_b());
            EXPECT_EQ(vector<double>({90, 90, 90}), eph.unc_theta());
            EXPECT_EQ(vector<double>({0, 1, 2}), eph.rh());
            EXPECT_EQ(vector<double>({1, 0, 1}), eph.delta());
            EXPECT_EQ(vector<double>({180, 0, 90}), eph.phase());
            EXPECT_EQ(vector<double>({0, 180, 80}), eph.selong());
            EXPECT_EQ(vector<double>({0, 30, 90}), eph.true_anomaly());
            EXPECT_EQ(vector<double>({0, 0, 0}), eph.sangle());
            EXPECT_EQ(vector<double>({10, 20, 30}), eph.vangle());
            EXPECT_EQ(vector<double>({-1, 5, 10}), eph.vmag());

            // initialize with invalid mjd order
            EXPECT_THROW(Ephemeris(encke, {data[1], data[0]}), std::runtime_error);
        }

        TEST_F(EphemerisTest, EphemerisVertex)
        {
            Ephemeris eph = Ephemeris(MovingTarget(), data);
            for (int i = 0; i < data.size(); i++)
                EXPECT_EQ(eph.vertex(i), eph.vertex(i - 3));

            EXPECT_THROW(eph.vertex(3), std::runtime_error);
            EXPECT_THROW(eph.vertex(-4), std::runtime_error);

            EXPECT_EQ(eph.vertex(0), data[0].as_s2point());
            EXPECT_EQ(eph.vertex(1), data[1].as_s2point());
            EXPECT_EQ(eph.vertex(2), data[2].as_s2point());
        }

        TEST_F(EphemerisTest, EphemerisStreamOperator)
        {
            std::stringstream s;
            Ephemeris eph(encke, data);
            s << eph;
            EXPECT_EQ(
                s.str(),
                "     mjd       tmtp        ra       dec      rh   delta    phase   selong  true_anomaly  sangle  vangle   unc_a  unc_b  unc_th    vmag\n"
                "--------  ---------  --------  --------  ------  ------  -------  -------  ------------  ------  ------  ------  -----  ------  ------\n"
                "0.000000  10.000000  1.000000  0.000000  0.0000  1.0000  180.000    0.000         0.000   0.000  10.000   1.000  0.100  90.000  -1.000\n"
                "1.000000  11.000000  2.000000  0.000000  1.0000  0.0000    0.000  180.000        30.000   0.000  20.000   5.000  0.500  90.000   5.000\n"
                "2.000000  12.000000  3.000000  0.000000  2.0000  1.0000   90.000   80.000        90.000   0.000  30.000  10.000  1.000  90.000  10.000\n");
        }

        TEST_F(EphemerisTest, EphemerisEquality)
        {
            Ephemeris eph(encke, data);
            Ephemeris same(encke, data);
            EXPECT_TRUE(eph == same);
            EXPECT_TRUE(same == eph);

            Ephemeris not_same = eph.segment(0);
            EXPECT_TRUE(eph != not_same);
            EXPECT_TRUE(not_same != eph);

            not_same = eph;
            not_same.target(MovingTarget("1P", 2));
            EXPECT_TRUE(eph != not_same);
            EXPECT_TRUE(not_same != eph);
        }

        TEST_F(EphemerisTest, EphemerisBracketOperator)
        {
            Ephemeris eph(encke, data);

            for (int i = 0; i < eph.num_vertices(); i++)
            {
                Ephemeris::Datum a = eph.data(i);
                EXPECT_EQ(a, eph[i].data(0));

                a = eph.data(i - eph.num_vertices());
                EXPECT_EQ(a, eph[i].data(0));
            }

            EXPECT_THROW(eph[eph.num_vertices()], std::runtime_error);
            EXPECT_THROW(eph[-eph.num_vertices() - 1], std::runtime_error);
        }

        TEST_F(EphemerisTest, EphemerisSlice)
        {
            Ephemeris eph(encke, data);
            Ephemeris subset(eph.slice(1));
            EXPECT_EQ(subset.num_vertices(), 2);
            EXPECT_EQ(subset.target(), encke);
            EXPECT_EQ(subset.data(0), eph.data(1));
            EXPECT_EQ(subset.data(1), eph.data(2));

            subset = eph.slice(0, 2);
            EXPECT_EQ(subset.num_vertices(), 2);
            EXPECT_EQ(subset.target(), encke);
            EXPECT_EQ(subset.data(0), eph.data(0));
            EXPECT_EQ(subset.data(1), eph.data(1));
        }

        TEST_F(EphemerisTest, EphemerisAppend)
        {
            Ephemeris eph(encke, data);
            Ephemeris a(encke, {{3, 13, 4, 1, 0, 0, 0, 3, 3, 45, 90, 50, 1, 2, 3}});
            eph.append(a);
            EXPECT_EQ(eph[3], a[0]);

            // append to an empty ephemeris
            Ephemeris b;
            EXPECT_THROW(b.append(a), std::runtime_error);
            b.target(a.target()); // need the object ids to match
            b.append(a);
            EXPECT_EQ(a, b);

            // don't append if the mjd is later
            EXPECT_THROW(a.append(eph), std::runtime_error);

            // don't append if the mjd is out of order
            EXPECT_THROW(a.append({{4}, {5}, {4.5}}), std::runtime_error);
        }

        TEST_F(EphemerisTest, EphemerisSegment)
        {
            Ephemeris eph(encke, data);

            // get single segment
            Ephemeris segment = eph.segment(1);
            Ephemeris expected(encke, {data[1], data[2]});
            EXPECT_EQ(segment, expected);

            // get all segments
            vector<Ephemeris> segments = eph.segments();
            for (int i = 0; i < eph.num_segments(); i++)
            {
                Ephemeris s = eph.segment(i);
                EXPECT_EQ(segments[i], s);

                s = eph.segment(i - eph.num_segments());
                EXPECT_EQ(segments[i], s);
            }

            // fail on invalid index
            EXPECT_THROW(eph.segment(eph.num_segments()), std::runtime_error);
            EXPECT_THROW(eph.segment(-eph.num_segments() - 1), std::runtime_error);
        }

        TEST_F(EphemerisTest, EphemerisAsPolyline)
        {
            Ephemeris eph(encke, data);
            S2Polyline polyline({data[0].as_s2point(), data[1].as_s2point(), data[2].as_s2point()});
            EXPECT_TRUE(eph.as_polyline().Equals(polyline));
        }

        TEST_F(EphemerisTest, EphemerisInterpolate)
        {
            Ephemeris eph(encke, data);

            Ephemeris interpolated = eph.interpolate(0.5);
            EXPECT_EQ(interpolated.target(), encke);
            EXPECT_EQ(interpolated.data(0).mjd, 0.5);
            EXPECT_EQ(interpolated.data(0).tmtp, 10.5);
            EXPECT_NEAR(interpolated.data(0).ra, 1.5, 1 * ARCSEC);
            EXPECT_NEAR(interpolated.data(0).dec, 0, 1 * ARCSEC);
            EXPECT_EQ(interpolated.data(0).unc_a, 3);
            EXPECT_NEAR(interpolated.data(0).unc_b, 0.3, 1e-8);
            EXPECT_EQ(interpolated.data(0).unc_theta, 90);
            EXPECT_EQ(interpolated.data(0).rh, 0.5);
            EXPECT_EQ(interpolated.data(0).delta, 0.5);
            EXPECT_EQ(interpolated.data(0).phase, 90);
            EXPECT_EQ(interpolated.data(0).selong, 90);
            EXPECT_EQ(interpolated.data(0).true_anomaly, 15);
            EXPECT_EQ(interpolated.data(0).sangle, 0);
            EXPECT_EQ(interpolated.data(0).vangle, 15);
            EXPECT_EQ(interpolated.data(0).vmag, 2);

            interpolated = eph.interpolate(1.5);
            EXPECT_NEAR(interpolated.data(0).ra, 2.5, 1 * ARCSEC);
            EXPECT_NEAR(interpolated.data(0).dec, 0, 1 * ARCSEC);

            // interpolate does not extrapolate
            EXPECT_THROW(eph.interpolate(-1), std::runtime_error);
            EXPECT_THROW(eph.interpolate(3), std::runtime_error);
        }

        TEST_F(EphemerisTest, EphemerisExtrapolate)
        {
            using Extrapolate = Ephemeris::Extrapolate;
            Ephemeris eph(encke, data);

            Ephemeris extrapolated = eph.extrapolate(5 * DEG, Extrapolate::BACKWARDS);
            EXPECT_EQ(extrapolated.target(), encke);
            EXPECT_NEAR(extrapolated.data(0).ra, -4, 1 * ARCSEC);
            EXPECT_NEAR(extrapolated.data(0).dec, 0, 1 * ARCSEC);

            extrapolated = eph.extrapolate(5 * DEG, Extrapolate::FORWARDS);
            EXPECT_NEAR(extrapolated.data(0).ra, 8, 1 * ARCSEC);
            EXPECT_NEAR(extrapolated.data(0).dec, 0, 1 * ARCSEC);
        }

        TEST_F(EphemerisTest, EphemerisSubsample)
        {
            Ephemeris eph(encke, data);

            Ephemeris subsample = eph.subsample(0.5, 0.75);
            EXPECT_EQ(subsample.target(), encke);
            EXPECT_EQ(subsample.num_segments(), 1);
            EXPECT_EQ(subsample[0], eph.interpolate(0.5));
            EXPECT_EQ(subsample[1], eph.interpolate(0.75));

            subsample = eph.subsample(0.5, 1.5);
            EXPECT_EQ(subsample.num_segments(), 2);
            EXPECT_EQ(subsample[0], eph.interpolate(0.5));
            EXPECT_EQ(subsample[1], eph[1]);
            EXPECT_EQ(subsample[2], eph.interpolate(1.5));

            subsample = eph.subsample(1, 2);
            EXPECT_EQ(subsample.num_segments(), 1);
            EXPECT_EQ(subsample[0], eph[1]);
            EXPECT_EQ(subsample[1], eph[2]);
        }

        void generate_expected_polygon(const S2LatLng &start, const S2LatLng &end, const double a, const double b, const double theta, S2Polygon &polygon)
        {
            // a, b, theta in radians

            // get the ellipses for the first and last points
            vector<S2LatLng> e0 = ellipse(16, start, a, b, theta);
            vector<S2LatLng> e1 = ellipse(16, end, a, b, theta);

            // our ephemeris varies by RA, and our padded region is elongated along Dec
            vector<S2LatLng> coords;
            // start with half of the first ellipse
            coords.insert(coords.end(), e0.begin() + 8, e0.end());
            coords.push_back(e0[0]);
            // append half the last ellipse
            coords.insert(coords.end(), e1.begin(), e1.begin() + 9);

            vector<S2Point> points(coords.size());
            std::transform(coords.begin(), coords.end(), points.begin(), [](S2LatLng c)
                           { return c.ToPoint(); });

            make_polygon(points, polygon);
        }

        TEST_F(EphemerisTest, EphemerisPad)
        {
            Ephemeris eph(encke, data);

            S2Polygon polygon;
            eph.pad(3600 * 2, 3600, 0, polygon);

            // which is equivalent to:
            S2Polygon polygon2;
            eph.pad(3600, 2 * 3600, polygon2);
            EXPECT_TRUE(polygon.BoundaryNear(polygon2, S1Angle::Radians(1 * ARCSEC)));

            S2Polygon expected;
            generate_expected_polygon(eph.data(0).as_s2latlng(), eph.data(2).as_s2latlng(), 2 * DEG, 1 * DEG, 0, expected);
            EXPECT_TRUE(polygon.BoundaryNear(expected, S1Angle::Radians(1 * ARCSEC)));

            EXPECT_THROW(eph.pad({1}, {1, 2, 3}, polygon), std::runtime_error);
            EXPECT_THROW(eph.pad({1, 2, 3}, {1, 2}, {0, 0, 0}, polygon), std::runtime_error);
            EXPECT_THROW(eph.pad(3600 * 90, 3600 * 90, polygon), std::runtime_error);
        }

        TEST_F(EphemerisTest, EphemerisAsPolygon)
        {
            Ephemeris eph(encke, {{0, 10, 1, 0, 10, 10, 90, 0, 1, 180, 0, 0, 0, 10, -1},
                                  {1, 11, 2, 0, 10, 10, 90, 1, 0, 0, 180, 30, 0, 20, 5},
                                  {2, 12, 3, 0, 10, 10, 90, 2, 1, 90, 80, 90, 0, 30, 10}});

            S2Polygon polygon;
            eph.as_polygon(polygon);

            S2Polygon expected;
            generate_expected_polygon(eph.data(0).as_s2latlng(), eph.data(2).as_s2latlng(), 0.1 * ARCSEC, 0.1 * ARCSEC, 0, expected);
            EXPECT_TRUE(polygon.BoundaryNear(expected, S1Angle::Radians(0.01 * ARCSEC)));

            eph.mutable_options()->use_uncertainty = true;
            eph.as_polygon(polygon);
            generate_expected_polygon(eph.data(0).as_s2latlng(), eph.data(2).as_s2latlng(), 10 * ARCSEC, 10 * ARCSEC, 0, expected);
            EXPECT_TRUE(polygon.BoundaryNear(expected, S1Angle::Radians(1 * ARCSEC)));
        }

        TEST_F(EphemerisTest, EphemerisAsJSON)
        {
            Ephemeris eph(encke, data);

            json::array vertices = eph.as_json();
            EXPECT_EQ(vertices.size(), 3);
            EXPECT_EQ(vertices.at(0).if_object()->at("mjd"), 0.);
            EXPECT_EQ(vertices.at(0).if_object()->at("tmtp"), 10.);
            EXPECT_EQ(vertices.at(0).if_object()->at("ra"), 1.);
            EXPECT_EQ(vertices.at(0).if_object()->at("dec"), 0.);
            EXPECT_EQ(vertices.at(0).if_object()->at("unc_a"), 1.);
            EXPECT_EQ(vertices.at(0).if_object()->at("unc_b"), 0.1);
            EXPECT_EQ(vertices.at(0).if_object()->at("unc_theta"), 90.);
            EXPECT_EQ(vertices.at(0).if_object()->at("rh"), 0.);
            EXPECT_EQ(vertices.at(0).if_object()->at("delta"), 1.);
            EXPECT_EQ(vertices.at(0).if_object()->at("phase"), 180.);
            EXPECT_EQ(vertices.at(0).if_object()->at("selong"), 0.);
            EXPECT_EQ(vertices.at(0).if_object()->at("true_anomaly"), 0.);
            EXPECT_EQ(vertices.at(0).if_object()->at("sangle"), 0.);
            EXPECT_EQ(vertices.at(0).if_object()->at("vangle"), 10.);
            EXPECT_EQ(vertices.at(0).if_object()->at("vmag"), -1.);
        }

    }
}