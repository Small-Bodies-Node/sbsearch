#include "config.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <optional>
#include <utility>
#include <set>
#include <string>
#include <sstream>
#include <s2/s2latlng.h>
#include <s2/s2metrics.h>
#include <s2/s2point.h>
#include <s2/s2polyline.h>
#include <s2/s2region_term_indexer.h>
#include <gtest/gtest.h>

#include "date.h"
#include "ephemeris.h"
#include "moving_target.h"
#include "query_info.h"
#include "sbsearch_testing.h"
#include "util/polygon.h"
#include "util/spherical.h"

using sbsearch::Ephemeris;
using std::optional;
using std::vector;

class EphemerisTest : public ::testing::Test
{
protected:
    Ephemeris::Data data{
        // mjd, tmtp, ra, dec, unc a, b, theta, rh, delta, phase, selong, true, sangle, vangle, vmag
        {0, 10, 1, 0, 100, 90, 1, 0.1, 90, 0, 1, 180, 0, 0, 0, 10, -1},
        {1, 11, 2, 0, 100, 90, 5, 0.5, 90, 1, 0, 0, 180, 30, 0, 20, 5},
        {2, 12, 3, 0, 90, 100, 10, 1.0, 90, 2, 1, 90, 80, 90, 0, 30, std::nullopt}};
    sbsearch::MovingTarget encke{"2P", 1, true};
};

namespace sbsearch::testing
{
    TEST(EphemerisTests, DatumEquality)
    {
        Ephemeris::Datum a{0, 10, 1, 0, 20, 90, 1, 0.1, 90, 0, 1, 180, 0, 0, 0, 10, -1};
        Ephemeris::Datum b{0, 10, 1, 0, 20, 90, 1, 0.1, 90, 0, 1, 180, 0, 0, 0, 10, -1};
        Ephemeris::Datum c{1, 11, 2, 0, 20, 100, 5, 0.5, 90, 1, 0, 0, 180, 30, 0, 20, 5};
        EXPECT_EQ(a, b);
        EXPECT_NE(a, c);
        EXPECT_NE(b, c);

        Ephemeris::Datum d{0, 10, 1, 0, 20, 90, 1, 0.1, 90, 0, 1, 180, 0, 0, 0, 10, -1};
        EXPECT_EQ(a, d);
        d.mjd = d.mjd.value() + 1;
        EXPECT_NE(a, d);

        d.mjd = a.mjd;
        EXPECT_EQ(a, d);
        d.tmtp = d.tmtp.value() + 1;
        EXPECT_NE(a, d);

        d.tmtp = a.tmtp;
        EXPECT_EQ(a, d);
        d.ra = d.ra.value() + 1;
        EXPECT_NE(a, d);

        d.ra = a.ra;
        EXPECT_EQ(a, d);
        d.dec = d.dec.value() + 1;
        EXPECT_NE(a, d);

        d.dec = a.dec;
        EXPECT_EQ(a, d);
        d.mu = d.mu.value() + 1;
        EXPECT_NE(a, d);

        d.mu = a.mu;
        EXPECT_EQ(a, d);
        d.mu_theta = d.mu_theta.value() + 1;
        EXPECT_NE(a, d);

        d.mu_theta = a.mu_theta;
        EXPECT_EQ(a, d);
        d.unc_a = d.unc_a.value() + 1;
        EXPECT_NE(a, d);

        d.unc_a = a.unc_a;
        EXPECT_EQ(a, d);
        d.unc_b = d.unc_b.value() + 1;
        EXPECT_NE(a, d);

        d.unc_b = a.unc_b;
        EXPECT_EQ(a, d);
        d.unc_theta = d.unc_theta.value() + 1;
        EXPECT_NE(a, d);

        d.unc_theta = a.unc_theta;
        EXPECT_EQ(a, d);
        d.rh = d.rh.value() + 1;
        EXPECT_NE(a, d);

        d.rh = a.rh;
        EXPECT_EQ(a, d);
        d.delta = d.delta.value() + 1;
        EXPECT_NE(a, d);

        d.delta = a.delta;
        EXPECT_EQ(a, d);
        d.phase = d.phase.value() + 1;
        EXPECT_NE(a, d);

        d.phase = a.phase;
        EXPECT_EQ(a, d);
        d.selong = d.selong.value() + 1;
        EXPECT_NE(a, d);

        d.selong = a.selong;
        EXPECT_EQ(a, d);
        d.true_anomaly = d.true_anomaly.value() + 1;
        EXPECT_NE(a, d);

        d.true_anomaly = a.true_anomaly;
        EXPECT_EQ(a, d);
        d.sangle = d.sangle.value() + 1;
        EXPECT_NE(a, d);

        d.sangle = a.sangle;
        EXPECT_EQ(a, d);
        d.vangle = d.vangle.value() + 1;
        EXPECT_NE(a, d);

        d.vangle = a.vangle;
        EXPECT_EQ(a, d);
        d.vmag = d.vmag.value() + 1;
        EXPECT_NE(a, d);
    }

    TEST(EphemerisTests, DatumRaDecFromPoint)
    {
        Ephemeris::Datum d;
        d.radec(S2LatLng::FromDegrees(0, 1).ToPoint());
        EXPECT_EQ(d.ra, 1);
        EXPECT_EQ(d.dec, 0);
    }

    TEST_F(EphemerisTest, DatumAsPoint)
    {
        EXPECT_EQ(data[0].as_s2point(), S2LatLng::FromDegrees(0, 1).ToPoint());
    }

    TEST_F(EphemerisTest, DatumAsJSON)
    {
        json::object obj = data[0].as_json();
        EXPECT_EQ(obj["mjd"], 0.);
        EXPECT_EQ(obj["tmtp"], 10.);
        EXPECT_EQ(obj["ra"], 1.);
        EXPECT_EQ(obj["dec"], 0.);
        EXPECT_EQ(obj["mu"], 100.);
        EXPECT_EQ(obj["mu_theta"], 90.);
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

    TEST_F(EphemerisTest, Init)
    {
        Ephemeris eph;
        EXPECT_EQ(eph.num_segments(), 0);
        EXPECT_EQ(eph.num_vertices(), 0);

        // single point ephemeris, number of segments should be 0
        eph = Ephemeris(encke, {data[0]});
        EXPECT_EQ(eph.num_segments(), 0);
        EXPECT_EQ(eph.num_vertices(), 1);

        eph = Ephemeris(encke, data);
        EXPECT_EQ(vector<optional<double>>({0, 1, 2}), eph.mjd());
        EXPECT_EQ(vector<optional<double>>({10, 11, 12}), eph.tmtp());
        EXPECT_EQ(vector<optional<double>>({1, 2, 3}), eph.ra());
        EXPECT_EQ(vector<optional<double>>({0, 0, 0}), eph.dec());
        EXPECT_EQ(vector<optional<double>>({100, 100, 90}), eph.mu());
        EXPECT_EQ(vector<optional<double>>({90, 90, 100}), eph.mu_theta());
        EXPECT_EQ(vector<optional<double>>({1, 5, 10}), eph.unc_a());
        EXPECT_EQ(vector<optional<double>>({0.1, 0.5, 1.0}), eph.unc_b());
        EXPECT_EQ(vector<optional<double>>({90, 90, 90}), eph.unc_theta());
        EXPECT_EQ(vector<optional<double>>({0, 1, 2}), eph.rh());
        EXPECT_EQ(vector<optional<double>>({1, 0, 1}), eph.delta());
        EXPECT_EQ(vector<optional<double>>({180, 0, 90}), eph.phase());
        EXPECT_EQ(vector<optional<double>>({0, 180, 80}), eph.selong());
        EXPECT_EQ(vector<optional<double>>({0, 30, 90}), eph.true_anomaly());
        EXPECT_EQ(vector<optional<double>>({0, 0, 0}), eph.sangle());
        EXPECT_EQ(vector<optional<double>>({10, 20, 30}), eph.vangle());
        EXPECT_EQ(vector<optional<double>>({-1, 5, std::nullopt}), eph.vmag());

        // initialize with invalid mjd order
        EXPECT_THROW(Ephemeris(encke, {data[1], data[0]}), std::runtime_error);
    }

    TEST_F(EphemerisTest, Vertex)
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

    TEST_F(EphemerisTest, StreamOperator)
    {
        std::stringstream s;
        Ephemeris eph(encke, data);
        s << eph;
        EXPECT_EQ(
            s.str(),
            "     mjd       tmtp        ra       dec      mu  mu_theta      rh   delta    phase   selong  true_anomaly  sangle  vangle   unc_a  unc_b  unc_th    vmag\n"
            "--------  ---------  --------  --------  ------  --------  ------  ------  -------  -------  ------------  ------  ------  ------  -----  ------  ------\n"
            "0.000000  10.000000  1.000000  0.000000  100.00    90.000  0.0000  1.0000  180.000    0.000         0.000   0.000  10.000   1.000  0.100  90.000  -1.000\n"
            "1.000000  11.000000  2.000000  0.000000  100.00    90.000  1.0000  0.0000    0.000  180.000        30.000   0.000  20.000   5.000  0.500  90.000   5.000\n"
            "2.000000  12.000000  3.000000  0.000000   90.00   100.000  2.0000  1.0000   90.000   80.000        90.000   0.000  30.000  10.000  1.000  90.000    null\n");

        eph.format.date = DateFormat::CALENDAR;
        s.str("");
        s << eph;
        EXPECT_EQ(
            s.str(),
            "               date       tmtp        ra       dec      mu  mu_theta      rh   delta    phase   selong  true_anomaly  sangle  vangle   unc_a  unc_b  unc_th    vmag\n"
            "-------------------  ---------  --------  --------  ------  --------  ------  ------  -------  -------  ------------  ------  ------  ------  -----  ------  ------\n"
            "1858-11-17 00:00:00  10.000000  1.000000  0.000000  100.00    90.000  0.0000  1.0000  180.000    0.000         0.000   0.000  10.000   1.000  0.100  90.000  -1.000\n"
            "1858-11-18 00:00:00  11.000000  2.000000  0.000000  100.00    90.000  1.0000  0.0000    0.000  180.000        30.000   0.000  20.000   5.000  0.500  90.000   5.000\n"
            "1858-11-19 00:00:00  12.000000  3.000000  0.000000   90.00   100.000  2.0000  1.0000   90.000   80.000        90.000   0.000  30.000  10.000  1.000  90.000    null\n");
    }

    TEST_F(EphemerisTest, Equality)
    {
        Ephemeris eph(encke, data);
        Ephemeris same(encke, data);
        EXPECT_TRUE(eph == same);
        EXPECT_TRUE(same == eph);

        Ephemeris not_same = eph.segment(0);
        EXPECT_TRUE(eph != not_same);
        EXPECT_TRUE(not_same != eph);

        not_same = eph;
        not_same.target(MovingTarget("1P", 2, true));
        EXPECT_TRUE(eph != not_same);
        EXPECT_TRUE(not_same != eph);
    }

    TEST_F(EphemerisTest, BracketOperator)
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

    TEST_F(EphemerisTest, Slice)
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

    TEST_F(EphemerisTest, Append)
    {
        Ephemeris eph(encke, data);
        Ephemeris a(encke, {{3, 13, 10, 65, 4, 1, 0, 0, 0, 3, 3, 45, 90, 50, 1, 2, 3}});
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

    TEST_F(EphemerisTest, Segment)
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

    TEST_F(EphemerisTest, Split)
    {
        Ephemeris eph(encke, data);
        auto segments = eph.split(3, 10);
        EXPECT_EQ(segments.size(), 1);
        EXPECT_EQ(segments[0].num_segments(), 2);

        segments = eph.split(0.9, 10);
        EXPECT_EQ(segments.size(), 2);
        EXPECT_EQ(segments[0].num_segments(), 1);
        EXPECT_EQ(segments[1].num_segments(), 1);

        segments = eph.split(3, 0.9);
        EXPECT_EQ(segments.size(), 2);
        EXPECT_EQ(segments[0].num_segments(), 1);
        EXPECT_EQ(segments[1].num_segments(), 1);
    }

    TEST_F(EphemerisTest, AsPolyline)
    {
        Ephemeris eph(encke, data);
        S2Polyline polyline({data[0].as_s2point(), data[1].as_s2point(), data[2].as_s2point()});
        EXPECT_TRUE(eph.as_polyline().Equals(polyline));
    }

    TEST_F(EphemerisTest, Interpolate)
    {
        Ephemeris eph(encke, data);

        Ephemeris interpolated = eph.interpolate(0.5);
        EXPECT_EQ(interpolated.target(), encke);
        EXPECT_EQ(interpolated.data(0).mjd.value(), 0.5);
        EXPECT_EQ(interpolated.data(0).tmtp.value(), 10.5);
        EXPECT_NEAR(interpolated.data(0).ra.value(), 1.5, 1 * ARCSEC);
        EXPECT_NEAR(interpolated.data(0).dec.value(), 0, 1 * ARCSEC);
        EXPECT_NEAR(interpolated.data(0).mu.value(), 100, 1e-6);
        EXPECT_NEAR(interpolated.data(0).mu_theta.value(), 90, 1e-6);
        EXPECT_EQ(interpolated.data(0).unc_a.value(), 3);
        EXPECT_NEAR(interpolated.data(0).unc_b.value(), 0.3, 1e-8);
        EXPECT_EQ(interpolated.data(0).unc_theta.value(), 90);
        EXPECT_EQ(interpolated.data(0).rh.value(), 0.5);
        EXPECT_EQ(interpolated.data(0).delta.value(), 0.5);
        EXPECT_EQ(interpolated.data(0).phase.value(), 90);
        EXPECT_EQ(interpolated.data(0).selong.value(), 90);
        EXPECT_EQ(interpolated.data(0).true_anomaly.value(), 15);
        EXPECT_EQ(interpolated.data(0).sangle.value(), 0);
        EXPECT_EQ(interpolated.data(0).vangle.value(), 15);
        EXPECT_EQ(interpolated.data(0).vmag.value(), 2);

        interpolated = eph.interpolate(1.5);
        EXPECT_NEAR(interpolated.data(0).ra.value(), 2.5, 1 * ARCSEC);
        EXPECT_NEAR(interpolated.data(0).dec.value(), 0, 1 * ARCSEC);

        // interpolate does not extrapolate
        EXPECT_THROW(eph.interpolate(-1), std::runtime_error);
        EXPECT_THROW(eph.interpolate(3), std::runtime_error);
    }

    TEST_F(EphemerisTest, Extrapolate)
    {
        using Extrapolate = Ephemeris::Extrapolate;
        Ephemeris eph(encke, data);

        Ephemeris extrapolated = eph.extrapolate(5 * DEG, Extrapolate::BACKWARDS);
        EXPECT_EQ(extrapolated.target(), encke);
        EXPECT_NEAR(extrapolated.data(0).ra.value(), -4, 1 * ARCSEC);
        EXPECT_NEAR(extrapolated.data(0).dec.value(), 0, 1 * ARCSEC);

        extrapolated = eph.extrapolate(5 * DEG, Extrapolate::FORWARDS);
        EXPECT_NEAR(extrapolated.data(0).ra.value(), 8, 1 * ARCSEC);
        EXPECT_NEAR(extrapolated.data(0).dec.value(), 0, 1 * ARCSEC);
    }

    TEST_F(EphemerisTest, Subsample)
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
        vector<S2LatLng> e0 = util::ellipse(16, start, a, b, theta);
        vector<S2LatLng> e1 = util::ellipse(16, end, a, b, theta);

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

        util::make_polygon(points, polygon);
    }

    TEST_F(EphemerisTest, AsPolygons)
    {
        Ephemeris eph(encke, {{0, 10.0, 0.1, 0, 360. / 1400, 90.0, 100, 10, 90, 0, 1, 180, 0, 0, 0, 10, -1},
                              {1, 11.0, 0.2, 0, 360. / 1400, 90.0, 100, 10, 90, 1, 0, 0, 180, 30, 0, 20, 5},
                              {2, 12.0, 0.3, 0, 360. / 1400, 90.0, 100, 10, 45, 2, 1, 90, 80, 90, 0, 30, 10},
                              {3, 13.0, 0.4, 0, 360. / 1400, 90.0, 100, 10, 0, 2, 1, 90, 80, 90, 0, 30, 10}});

        eph.mutable_options()->use_uncertainty = true;
        auto polygons = eph.as_polygons();

        // // debugging
        // QueryInfo info;
        // for (auto &polygon : polygons)
        // {
        //     array<array<double, 2>, 4> vertices;
        //     for (int i = 0; i < 4; i++)
        //         vertices[i] = {S2LatLng::Longitude(polygon->loop(0)->vertex(i)).degrees(),
        //                        S2LatLng::Latitude(polygon->loop(0)->vertex(i)).degrees()};
        //     info.query_polygons.emplace_back(vertices);
        // }

        // vector<array<double, 5>> vector;
        // for (auto const &e : eph.data())
        //     vector.emplace_back(
        //         array<double, 5>{e.ra.value_or(1e99),
        //                          e.dec.value_or(1e99),
        //                          e.unc_a.value_or(1e99),
        //                          e.unc_b.value_or(1e99),
        //                          e.unc_theta.value_or(1e99)});

        // info.ephemeris_segments.emplace_back(std::move(vector));

        // std::ofstream os("test.json");
        // os << info;

        auto contains = [&polygons](S2Point p)
        {
            bool contained = false;
            for (auto const &polygon : polygons)
            {
                if (polygon->Contains(p))
                {
                    contained = true;
                    break;
                }
            }
            EXPECT_TRUE(contained);
        };

        // The polygons must contain all the ephemeris points and given
        // error ellipses.
        for (auto const &e : eph.data())
        {
            auto center = S2LatLng::FromDegrees(e.dec.value(), e.ra.value());
            contains(center.ToPoint());

            double a = e.unc_a.value() * ARCSEC;
            double b = e.unc_b.value() * ARCSEC;
            double theta = e.unc_theta.value() * DEG;
            for (auto const &point : util::ellipse(36, center, a, b, theta))
                contains(point.ToPoint());
        }

        // Sample the ephemeris.  The ephemeris must be contained.
        for (double mjd : {0.1, 0.5, 0.9, 1.3, 1.5, 1.8})
        {
            auto sample = eph.interpolate(mjd);
            contains(sample.vertex(0));
        }
    }

    TEST_F(EphemerisTest, AsJSON)
    {
        Ephemeris eph(encke, data);

        json::array vertices = eph.as_json();
        EXPECT_EQ(vertices.size(), 3);
        EXPECT_EQ(vertices.at(0).at("mjd"), 0.);
        EXPECT_EQ(vertices.at(0).at("tmtp"), 10.);
        EXPECT_EQ(vertices.at(0).at("ra"), 1.);
        EXPECT_EQ(vertices.at(0).at("dec"), 0.);
        EXPECT_EQ(vertices.at(0).at("mu"), 100.);
        EXPECT_EQ(vertices.at(0).at("mu_theta"), 90.);
        EXPECT_EQ(vertices.at(0).at("unc_a"), 1.);
        EXPECT_EQ(vertices.at(0).at("unc_b"), 0.1);
        EXPECT_EQ(vertices.at(0).at("unc_theta"), 90.);
        EXPECT_EQ(vertices.at(0).at("rh"), 0.);
        EXPECT_EQ(vertices.at(0).at("delta"), 1.);
        EXPECT_EQ(vertices.at(0).at("phase"), 180.);
        EXPECT_EQ(vertices.at(0).at("selong"), 0.);
        EXPECT_EQ(vertices.at(0).at("true_anomaly"), 0.);
        EXPECT_EQ(vertices.at(0).at("sangle"), 0.);
        EXPECT_EQ(vertices.at(0).at("vangle"), 10.);
        EXPECT_TRUE(vertices.at(2).at("vmag").is_null());
    }

}
