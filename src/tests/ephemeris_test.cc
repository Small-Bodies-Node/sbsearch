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

#include "config.h"
#include "constants.h"
#include "date.h"
#include "json.h"
#include "moving_target.h"
#include "polyline.h"
#include "polygons.h"
#include "query_info.h"
#include "sbsearch_testing.h"
#include "vertices.h"
#include "ephemeris/ephemeris.h"
#include "ephemeris/extrapolate.h"
#include "ephemeris/interpolate.h"
#include "ephemeris/parallax_offset.h"
#include "ephemeris/safe_sampler.h"
#include "ephemeris/split.h"
#include "ephemeris/subsample.h"
#include "util/optional.h"
#include "util/spherical.h"

using namespace sbsearch::ephemeris;

using std::optional;
using std::vector;

class EphemerisTest : public ::testing::Test
{
protected:
    Ephemeris::Data data{
        // ra = 1 + 0.1 * mjd + 0.02 *mjd**2 + 0.0005 * mjd**3
        // dec = -1.1 * mjd - 0.2 *mjd**2
        // mjd, tmtp, ra, dec, mu, mu_theta, unc a, b, theta, rh, delta, phase, selong, true, sangle, vangle, vmag
        {0, 10, 1.0000, 0.0, 100, 80, 1, 0.1, 90, 0, 1, 180, 0, 0, 0, 10, -1},
        {1, 11, 1.1205, -1.3, 90, 90, 5, 0.5, 90, 1, 0, 0, 180, 30, 0, 20, 5},
        {2, 12, 1.2840, -3.0, 80, 100, 10, 1.0, 90, 2, 1, 90, 80, 90, 0, 30, std::nullopt}};
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
        d.mjd = d.mjd + 1;
        EXPECT_NE(a, d);

        d.mjd = a.mjd;
        EXPECT_EQ(a, d);
        d.tmtp = d.tmtp.value() + 1;
        EXPECT_NE(a, d);

        d.tmtp = a.tmtp;
        EXPECT_EQ(a, d);
        d.ra = d.ra + 1;
        EXPECT_NE(a, d);

        d.ra = a.ra;
        EXPECT_EQ(a, d);
        d.dec = d.dec + 1;
        EXPECT_NE(a, d);

        d.dec = a.dec;
        EXPECT_EQ(a, d);
        d.mu = d.mu + 1;
        EXPECT_NE(a, d);

        d.mu = a.mu;
        EXPECT_EQ(a, d);
        d.mu_theta = d.mu_theta + 1;
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
        d.rh = d.rh + 1;
        EXPECT_NE(a, d);

        d.rh = a.rh;
        EXPECT_EQ(a, d);
        d.delta = d.delta + 1;
        EXPECT_NE(a, d);

        d.delta = a.delta;
        EXPECT_EQ(a, d);
        d.phase = d.phase + 1;
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
        EXPECT_EQ(make_point(data[0]), S2LatLng::FromDegrees(0, 1).ToPoint());
    }

    TEST_F(EphemerisTest, DatumAsJSON)
    {
        json::object obj(json::value_from(data[0]).as_object());
        EXPECT_EQ(obj["date"], 0.);
        EXPECT_EQ(obj["tmtp"], 10.);
        EXPECT_EQ(obj["ra"], 1.);
        EXPECT_EQ(obj["dec"], 0.);
        EXPECT_EQ(obj["mu"], 100.);
        EXPECT_EQ(obj["mu_theta"], 80.);
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
        EXPECT_EQ(vector<double>({0, 1, 2}), eph.mjd());
        EXPECT_EQ(vector<optional<double>>({10, 11, 12}), eph.tmtp());
        EXPECT_EQ(vector<double>({1, 1.1205, 1.2840}), eph.ra());
        EXPECT_EQ(vector<double>({0, -1.3, -3.0}), eph.dec());
        EXPECT_EQ(vector<double>({100, 90, 80}), eph.mu());
        EXPECT_EQ(vector<double>({80, 90, 100}), eph.mu_theta());
        EXPECT_EQ(vector<optional<double>>({1, 5, 10}), eph.unc_a());
        EXPECT_EQ(vector<optional<double>>({0.1, 0.5, 1.0}), eph.unc_b());
        EXPECT_EQ(vector<optional<double>>({90, 90, 90}), eph.unc_theta());
        EXPECT_EQ(vector<double>({0, 1, 2}), eph.rh());
        EXPECT_EQ(vector<double>({1, 0, 1}), eph.delta());
        EXPECT_EQ(vector<double>({180, 0, 90}), eph.phase());
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

        EXPECT_EQ(eph.vertex(0), make_point(data[0]));
        EXPECT_EQ(eph.vertex(1), make_point(data[1]));
        EXPECT_EQ(eph.vertex(2), make_point(data[2]));
    }

    TEST_F(EphemerisTest, StreamOperator)
    {
        std::stringstream s;
        Ephemeris eph(encke, data);
        s << eph;
        EXPECT_EQ(
            s.str(),
            "    date       tmtp        ra        dec      mu  mu_theta      rh   delta    phase   selong  true_anomaly  sangle  vangle   unc_a  unc_b  unc_th    vmag\n"
            "--------  ---------  --------  ---------  ------  --------  ------  ------  -------  -------  ------------  ------  ------  ------  -----  ------  ------\n"
            "0.000000  10.000000  1.000000   0.000000  100.00    80.000  0.0000  1.0000  180.000    0.000         0.000   0.000  10.000   1.000  0.100  90.000  -1.000\n"
            "1.000000  11.000000  1.120500  -1.300000   90.00    90.000  1.0000  0.0000    0.000  180.000        30.000   0.000  20.000   5.000  0.500  90.000   5.000\n"
            "2.000000  12.000000  1.284000  -3.000000   80.00   100.000  2.0000  1.0000   90.000   80.000        90.000   0.000  30.000  10.000  1.000  90.000    null\n");

        eph.format.date = Date::Format::Calendar;
        s.str("");
        s << eph;
        EXPECT_EQ(
            s.str(),
            "               date       tmtp        ra        dec      mu  mu_theta      rh   delta    phase   selong  true_anomaly  sangle  vangle   unc_a  unc_b  unc_th    vmag\n"
            "-------------------  ---------  --------  ---------  ------  --------  ------  ------  -------  -------  ------------  ------  ------  ------  -----  ------  ------\n"
            "1858-11-17 00:00:00  10.000000  1.000000   0.000000  100.00    80.000  0.0000  1.0000  180.000    0.000         0.000   0.000  10.000   1.000  0.100  90.000  -1.000\n"
            "1858-11-18 00:00:00  11.000000  1.120500  -1.300000   90.00    90.000  1.0000  0.0000    0.000  180.000        30.000   0.000  20.000   5.000  0.500  90.000   5.000\n"
            "1858-11-19 00:00:00  12.000000  1.284000  -3.000000   80.00   100.000  2.0000  1.0000   90.000   80.000        90.000   0.000  30.000  10.000  1.000  90.000    null\n");
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
        Ephemeris eph(encke, data, {Date::Format::Calendar});
        Ephemeris subset(eph.slice(1));
        EXPECT_EQ(subset.num_vertices(), 2);
        EXPECT_EQ(subset.target(), encke);
        EXPECT_EQ(subset.data(0), eph.data(1));
        EXPECT_EQ(subset.data(1), eph.data(2));
        EXPECT_EQ(subset.format.date, Date::Format::Calendar);

        eph.format.date = Date::Format::MJD;
        subset = eph.slice(0, 2);
        EXPECT_EQ(subset.num_vertices(), 2);
        EXPECT_EQ(subset.target(), encke);
        EXPECT_EQ(subset.data(0), eph.data(0));
        EXPECT_EQ(subset.data(1), eph.data(1));
        EXPECT_EQ(subset.format.date, Date::Format::MJD);
    }

    TEST_F(EphemerisTest, Append)
    {
        Ephemeris eph(encke, data);
        Ephemeris a(encke, {{3, 13, 10, 65, 4, 1, 0, 0, 0, 3, 3, 45, 90, 50, 1, 2, 3}});
        eph.append(a);
        EXPECT_EQ(eph[3], a[0]);
        eph.append(std::move(a));
        EXPECT_EQ(eph[4], a[0]);

        // append to an empty ephemeris
        Ephemeris b;
        EXPECT_THROW(b.append(a), std::runtime_error);
        b.target(a.target()); // need the object ids to match
        b.append(a);
        EXPECT_EQ(a, b);

        // don't append if the mjd is later
        EXPECT_THROW(a.append(eph), std::runtime_error);

        // don't append if the mjd is out of order
        EXPECT_THROW(a.append(Ephemeris::Data{{4}, {5}, {4.5}}), std::runtime_error);
    }

    TEST_F(EphemerisTest, Segment)
    {
        Ephemeris eph(encke, data);

        // get single segment
        Ephemeris segment = eph.segment(1);
        Ephemeris expected(encke, {data[1], data[2]});
        EXPECT_EQ(segment, expected);

        // test date format
        EXPECT_EQ(segment.format.date, Date::Format::MJD);
        eph.format.date = Date::Format::Calendar;
        segment = eph.segment(1);
        EXPECT_EQ(segment.format.date, Date::Format::Calendar);

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
        auto segments = split(eph, 3, 10);
        EXPECT_EQ(segments.size(), 1);
        EXPECT_EQ(segments[0].num_segments(), 2);

        eph.format.date = Date::Format::Calendar;
        segments = split(eph, 0.9, 10);
        EXPECT_EQ(segments.size(), 2);
        EXPECT_EQ(segments[0].num_segments(), 1);
        EXPECT_EQ(segments[1].num_segments(), 1);
        EXPECT_EQ(segments[0].format.date, Date::Format::Calendar);

        segments = split(eph, 3, 0.9);
        EXPECT_EQ(segments.size(), 2);
        EXPECT_EQ(segments[0].num_segments(), 1);
        EXPECT_EQ(segments[1].num_segments(), 1);
    }

    TEST_F(EphemerisTest, AsPolyline)
    {
        Ephemeris eph(encke, data);
        S2Polyline polyline({make_point(data[0]), make_point(data[1]), make_point(data[2])});
        EXPECT_TRUE(make_polyline(eph).Equals(polyline));
    }

    TEST_F(EphemerisTest, Interpolation)
    {
        Ephemeris eph(encke, data);

        Ephemeris::Datum interpolated = interpolate(eph.data(), 0.5);
        EXPECT_EQ(interpolated.mjd, 0.5);
        EXPECT_EQ(interpolated.tmtp.value(), 10.5);
        EXPECT_NEAR(interpolated.ra * DEG, 1.0550625 * DEG, 1 * ARCSEC);
        EXPECT_NEAR(interpolated.dec * DEG, -0.6 * DEG, 1 * ARCSEC);
        EXPECT_NEAR(interpolated.mu, 95, 1e-6);
        EXPECT_NEAR(interpolated.mu_theta, 85, 1e-6);
        EXPECT_EQ(interpolated.unc_a.value(), 2.875);
        EXPECT_NEAR(interpolated.unc_b.value(), 0.2875, 1e-8);
        EXPECT_EQ(interpolated.unc_theta.value(), 90);
        EXPECT_EQ(interpolated.rh, 0.5);
        EXPECT_EQ(interpolated.delta, 0.25);
        EXPECT_EQ(interpolated.phase, 56.25);
        EXPECT_EQ(interpolated.selong.value(), 125);
        EXPECT_EQ(interpolated.true_anomaly.value(), 11.25);
        EXPECT_EQ(interpolated.sangle.value(), 0);
        EXPECT_EQ(interpolated.vangle, 15);
        EXPECT_FALSE(interpolated.vmag.has_value());

        interpolated = interpolate(eph.data(), 1.5);
        EXPECT_NEAR(interpolated.ra * DEG, 1.1966875 * DEG, 1 * ARCSEC);
        EXPECT_NEAR(interpolated.dec * DEG, -2.1 * DEG, 1 * ARCSEC);

        // test target and format
        Ephemeris interpolated_eph = interpolate(eph, 1.5);
        EXPECT_EQ(interpolated_eph.target(), encke);
        EXPECT_EQ(interpolated_eph.format.date, Date::Format::MJD);

        eph.format.date = Date::Format::Calendar;
        interpolated_eph = interpolate(eph, 1.5);
        EXPECT_EQ(interpolated_eph.format.date, Date::Format::Calendar);

        // interpolate does not extrapolate
        EXPECT_THROW(interpolate(eph, -1), std::runtime_error);
        EXPECT_THROW(interpolate(eph, 3), std::runtime_error);
    }

    TEST_F(EphemerisTest, Extrapolate)
    {
        Ephemeris eph(encke, {{0, 10, 1, 0.0, 100, 80, 1, 0.1, 90, 0, 1, 180, 0, 0, 0, 10, -1},
                              {1, 11, 2, 0.0, 90, 90, 5, 0.5, 90, 1, 0, 0, 180, 30, 0, 20, 5}});

        Ephemeris::Datum extrapolated = extrapolate(eph.data(), 5 * DEG, Extrapolate::BACKWARDS);
        EXPECT_NEAR(extrapolated.ra, -4, 1 * ARCSEC / DEG);
        EXPECT_NEAR(extrapolated.dec, 0, 1 * ARCSEC / DEG);

        extrapolated = extrapolate(eph.data(), 5 * DEG, Extrapolate::FORWARDS);
        EXPECT_NEAR(extrapolated.ra, 7, 1 * ARCSEC / DEG);
        EXPECT_NEAR(extrapolated.dec, 0, 1 * ARCSEC / DEG);

        // test target and format
        Ephemeris extrapolated_eph = extrapolate(eph, 5 * DEG, Extrapolate::BACKWARDS);
        EXPECT_EQ(extrapolated_eph.target(), encke);
        EXPECT_EQ(extrapolated_eph.format.date, Date::Format::MJD);

        eph.format.date = Date::Format::Calendar;
        extrapolated_eph = extrapolate(eph, 5 * DEG, Extrapolate::BACKWARDS);
        EXPECT_EQ(extrapolated_eph.format.date, Date::Format::Calendar);
    }

    TEST_F(EphemerisTest, Subsample)
    {
        Ephemeris eph(encke, data);

        Ephemeris subsampled = subsample(eph, 0.5, 0.75);
        EXPECT_EQ(subsampled.target(), encke);
        EXPECT_EQ(subsampled.num_segments(), 1);
        EXPECT_EQ(subsampled.data(0), interpolate(eph.data(), 0.5));
        EXPECT_EQ(subsampled.data(1), interpolate(eph.data(), 0.75));
        EXPECT_EQ(subsampled.format.date, Date::Format::MJD);

        eph.format.date = Date::Format::Calendar;
        subsampled = subsample(eph, 0.5, 1.5);
        EXPECT_EQ(subsampled.num_segments(), 2);
        EXPECT_EQ(subsampled.data(0), interpolate(eph.data(), 0.5));
        EXPECT_EQ(subsampled[1], eph[1]);
        EXPECT_EQ(subsampled.data(2), interpolate(eph.data(), 1.5));
        EXPECT_EQ(subsampled.format.date, Date::Format::Calendar);

        subsampled = subsample(eph, 1, 2);
        EXPECT_EQ(subsampled.num_segments(), 1);
        EXPECT_EQ(subsampled[0], eph[1]);
        EXPECT_EQ(subsampled[1], eph[2]);

        // subsample empty ephemeris
        eph = Ephemeris();
        subsampled = subsample(eph, 0, 1);
        EXPECT_EQ(subsampled.num_vertices(), 0);
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

        make_polygon(points, polygon);
    }

    TEST_F(EphemerisTest, AsPolygons)
    {
        Ephemeris eph(encke, {{0, 10.0, 0.1, 0, 360. / 1400, 90.0, 100, 10, 90, 0, 1, 180, 0, 0, 0, 10, -1},
                              {1, 11.0, 0.2, 0, 360. / 1400, 90.0, 100, 10, 90, 1, 0, 0, 180, 30, 0, 20, 5},
                              {2, 12.0, 0.3, 0, 360. / 1400, 90.0, 100, 10, 45, 2, 1, 90, 80, 90, 0, 30, 10},
                              {3, 13.0, 0.4, 0, 360. / 1400, 90.0, 100, 10, 0, 2, 1, 90, 80, 90, 0, 30, 10}});

        auto polygons = make_polygons(eph, true, 0);

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
            auto center = S2LatLng::FromDegrees(e.dec, e.ra);
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
            auto sample = interpolate(eph.data(), mjd);
            contains(make_point(sample));
        }
    }

    TEST_F(EphemerisTest, AsJSON)
    {
        Ephemeris eph(encke, data);

        json::object obj = json::value_from(eph).as_object();
        EXPECT_EQ(obj["target"], json::value_from(eph.target()).as_object());

        json::array vertices = obj["data"].as_array();
        EXPECT_EQ(vertices.size(), 3);
        EXPECT_EQ(vertices.at(0).at("date"), 0.);
        EXPECT_EQ(vertices.at(0).at("tmtp"), 10.);
        EXPECT_EQ(vertices.at(0).at("ra"), 1.);
        EXPECT_EQ(vertices.at(0).at("dec"), 0.);
        EXPECT_EQ(vertices.at(0).at("mu"), 100.);
        EXPECT_EQ(vertices.at(0).at("mu_theta"), 80.);
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

        eph.format.date = Date::Format::Calendar;
        obj = json::value_from(eph).as_object();
        vertices = obj["data"].as_array();
        EXPECT_EQ(vertices.at(0).at("date"), "1858-11-17 00:00:00");
    }

    TEST_F(EphemerisTest, SafeSampler)
    {
        Ephemeris eph(encke, data);
        SafeSampler sampler;
        sampler.append(eph);
        Ephemeris result = sampler.subsample(0.5, 0.6);
        EXPECT_EQ(result.target(), eph.target());
    }
}
