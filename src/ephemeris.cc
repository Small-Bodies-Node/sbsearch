#include "config.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <stdexcept>
#include <tuple>
#include <vector>
#include <s2/s1angle.h>
#include <s2/s1chord_angle.h>
#include <s2/s2convex_hull_query.h>
#include <s2/s2edge_distances.h>
#include <s2/s2latlng.h>
#include <s2/s2point.h>
#include <s2/s2polyline.h>
#include <s2/s2region_term_indexer.h>

#include "ephemeris.h"
#include "exceptions.h"
#include "observatory.h"
#include "table.h"
#include "util/math.h"
#include "util/spherical.h"

using sbsearch::table::Table;
using std::cerr;
using std::cout;
using std::endl;
using std::unique_ptr;
using std::vector;

namespace sbsearch
{
    bool Ephemeris::Datum::operator==(const Ephemeris::Datum &other) const
    {
        return (std::tie(mjd, tmtp, ra, dec, unc_a, unc_b, unc_theta,
                         rh, delta, phase, selong, true_anomaly,
                         sangle, vangle, vmag) ==
                std::tie(other.mjd, other.tmtp, other.ra, other.dec, other.unc_a, other.unc_b, other.unc_theta,
                         other.rh, other.delta, other.phase, other.selong, other.true_anomaly,
                         other.sangle, other.vangle, other.vmag));
    }

    bool Ephemeris::Datum::operator!=(const Ephemeris::Datum &other) const
    {
        return !(*this == other);
    }

    json::object Ephemeris::Datum::as_json()
    {
        json::object datum;
        if (mjd)
            datum["mjd"] = mjd.value();

        if (tmtp)
            datum["tmtp"] = tmtp.value();

        if (ra)
            datum["ra"] = ra.value();

        if (dec)
            datum["dec"] = dec.value();

        if (unc_a)
            datum["unc_a"] = unc_a.value();

        if (unc_b)
            datum["unc_b"] = unc_b.value();

        if (unc_theta)
            datum["unc_theta"] = unc_theta.value();

        if (rh)
            datum["rh"] = rh.value();

        if (delta)
            datum["delta"] = delta.value();

        if (phase)
            datum["phase"] = phase.value();

        if (selong)
            datum["selong"] = selong.value();

        if (true_anomaly)
            datum["true_anomaly"] = true_anomaly.value();

        if (sangle)
            datum["sangle"] = sangle.value();

        if (vangle)
            datum["vangle"] = vangle.value();

        if (vmag)
            datum["vmag"] = vmag.value();

        return datum;
    }

    Ephemeris::Ephemeris(const MovingTarget target, Data data)
    {
        num_vertices_ = data.size();
        num_segments_ = (num_vertices_ == 0) ? 0 : (num_vertices_ - 1);

        target_ = target;
        data_ = Data(data);

        isValid();
    }

    const Ephemeris Ephemeris::operator[](const int k) const
    {
        Ephemeris eph = Ephemeris(target_, {data(k)});
        return eph;
    }

    const Ephemeris Ephemeris::slice(const int start) const
    {
        const int i = normalize_index(start, num_vertices_);
        Data subset(data_.begin() + i, data_.end());
        return Ephemeris(target_, subset);
    }

    const Ephemeris Ephemeris::slice(const int start, const int stop) const
    {
        const int i = normalize_index(start, num_vertices_);
        const int j = normalize_index(stop, num_vertices_ + 1);

        if (i > j)
            throw EphemerisError("start cannot be greater than stop.");

        Data subset(data_.begin() + i, data_.begin() + j);
        return Ephemeris(target_, subset);
    }

    bool Ephemeris::isValid() const
    {
        if (!util::is_increasing(mjd()))
            throw std::runtime_error("mjd must be monotonically increasing.");

        return true;
    }

    std::ostream &operator<<(std::ostream &os, const Ephemeris &ephemeris)
    {
        Table table;
        table.add_column("mjd", "%.6lf", ephemeris.mjd());
        table.add_column("tmtp", "%.6lf", ephemeris.tmtp());
        table.add_column("ra", "%.6lf", ephemeris.ra());
        table.add_column("dec", "%.6lf", ephemeris.dec());
        table.add_column("rh", "%.4f", ephemeris.rh());
        table.add_column("delta", "%.4f", ephemeris.delta());
        table.add_column("phase", "%.3f", ephemeris.phase());
        table.add_column("selong", "%.3f", ephemeris.selong());
        table.add_column("true_anomaly", "%.3f", ephemeris.true_anomaly());
        table.add_column("sangle", "%.3f", ephemeris.sangle());
        table.add_column("vangle", "%.3f", ephemeris.vangle());
        table.add_column("unc_a", "%.3f", ephemeris.unc_a());
        table.add_column("unc_b", "%.3f", ephemeris.unc_b());
        table.add_column("unc_th", "%.3f", ephemeris.unc_theta());
        table.add_column("vmag", "%.3f", ephemeris.vmag());

        os << table;
        return os;
    }

    bool
    Ephemeris::operator==(const Ephemeris &other) const
    {
        return ((target_ == other.target()) & (data_ == other.data()));
    }

    int Ephemeris::num_vertices() const
    {
        return num_vertices_;
    }

    const Ephemeris::Datum &Ephemeris::data(const int k) const
    {
        const int i = normalize_index(k, num_vertices_);
        return data_[i];
    }

    vector<double> Ephemeris::mjd() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.mjd.value_or(-1); });
        return result;
    }

    vector<double> Ephemeris::tmtp() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.tmtp.value_or(-1); });
        return result;
    }

    vector<double> Ephemeris::ra() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.ra.value_or(-1); });
        return result;
    }

    vector<double> Ephemeris::dec() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.dec.value_or(-1); });
        return result;
    }

    vector<double> Ephemeris::unc_a() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.unc_a.value_or(-1); });
        return result;
    }

    vector<double> Ephemeris::unc_b() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.unc_b.value_or(-1); });
        return result;
    }

    vector<double> Ephemeris::unc_theta() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.unc_theta.value_or(0); });
        return result;
    }

    vector<double> Ephemeris::rh() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.rh.value_or(-1); });
        return result;
    }

    vector<double> Ephemeris::delta() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.delta.value_or(-1); });
        return result;
    }

    vector<double> Ephemeris::phase() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.phase.value_or(-1); });
        return result;
    }

    vector<double> Ephemeris::selong() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.selong.value_or(-1); });
        return result;
    }

    vector<double> Ephemeris::true_anomaly() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.true_anomaly.value_or(-1); });
        return result;
    }

    vector<double> Ephemeris::sangle() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.sangle.value_or(-1); });
        return result;
    }

    vector<double> Ephemeris::vangle() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.vangle.value_or(-1); });
        return result;
    }

    vector<double> Ephemeris::vmag() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.vmag.value_or(99); });
        return result;
    }

    S2Point Ephemeris::vertex(const int k) const
    {
        return data(k).as_s2point();
    }

    vector<S2Point> Ephemeris::vertices() const
    {
        vector<S2Point> result(num_vertices_);
        for (int i = 0; i < num_vertices_; i++)
            result[i] = vertex(i);
        return result;
    }

    int Ephemeris::num_segments() const
    {
        return num_segments_;
    }

    void Ephemeris::append(const Data &new_data)
    {
        if (num_vertices_ != 0)
            if (data(-1).mjd > new_data[0].mjd)
                throw std::runtime_error("Attempting to append an ephemeris with an earlier mjd.");

        // check that new_data's time axis is OK
        auto i = std::adjacent_find(new_data.begin(), new_data.end(),
                                    [](const Datum &a, const Datum &b)
                                    { return a.mjd > b.mjd; });
        if (i != new_data.end())
            throw std::runtime_error("mjd must be monotonically increasing.");

        data_.insert(data_.end(), new_data.begin(), new_data.end());
        num_vertices_ = data_.size();
        num_segments_ = (num_vertices_ == 0) ? 0 : (num_vertices_ - 1);
    }

    void Ephemeris::append(const Ephemeris &eph)
    {
        if (eph.target().moving_target_id() != target_.moving_target_id())
            throw std::runtime_error("Attempted to append an ephemeris with a different object ID.");

        append(eph.data());
    }

    Ephemeris Ephemeris::segment(const int k) const
    {
        const int i = normalize_index(k, num_segments_);
        return Ephemeris(target_, {data_[i], data_[i + 1]});
    }

    vector<Ephemeris> Ephemeris::segments() const
    {
        vector<Ephemeris> eph(num_segments_);
        for (int i = 0; i < num_segments_; i++)
            eph[i] = segment(i);
        return eph;
    }

    vector<Ephemeris> Ephemeris::split(double length, double time) const
    {
        if (num_vertices_ <= 1)
            return {};

        vector<Ephemeris> segments;
        segments.reserve(std::ceil(as_polyline().GetLength().degrees() / length));
        double arc = 0, period = 0;
        int start = 0;
        for (int i = 0; i < num_segments_; i++)
        {
            Ephemeris segment_ = segment(i);
            arc += segment_.as_polyline().GetLength().degrees();
            period += segment_.data(1).mjd.value() - segment_.data(0).mjd.value();

            if ((arc >= length) | (period >= time) | (i == num_segments_ - 1))
            {
                segments.push_back(slice(start, i + 2));
                arc = 0;
                period = 0;
                start = i + 1;
            }
        }
        return segments;
    }

    S2Polyline Ephemeris::as_polyline() const
    {
        return S2Polyline(vertices());
    }

    Ephemeris Ephemeris::parallax_offset(const Observatory &observatory)
    {
        Data new_data(data_);
        for (int i = 0; i < data_.size(); i++)
        {
            const S2LatLng coords = observatory.parallax(
                new_data[i].as_s2latlng(),
                new_data[i].mjd.value(),
                new_data[i].delta.value());
            new_data[i].radec(coords.Normalized().ToPoint());
        }
        return Ephemeris(target_, new_data);
    }

    Ephemeris Ephemeris::interpolate(const double mjd) const
    {
        if ((mjd < data_.front().mjd) | (mjd > data_.back().mjd))
            throw std::runtime_error("Interpolation beyond ephemeris time range: ");

        // find the nearest segment
        auto end = std::find_if(data_.begin(), data_.end(), [mjd](Datum d)
                                { return (d.mjd.value() > mjd); });
        auto start = end - 1;

        // length of segment in time
        double dt = (*end).mjd.value() - (*start).mjd.value();

        // interpolate to this fraction
        double frac = (mjd - (*start).mjd.value()) / dt;

        // this is the line we will interpolate
        S2Polyline segment(vector<S2Point>({(*start).as_s2point(), (*end).as_s2point()}));

        Datum d;
        d.mjd = util::interp((*start).mjd, (*end).mjd, frac);
        d.tmtp = util::interp((*start).tmtp, (*end).tmtp, frac);
        d.radec(segment.Interpolate(frac));
        d.unc_a = util::interp((*start).unc_a, (*end).unc_a, frac);
        d.unc_b = util::interp((*start).unc_b, (*end).unc_b, frac);
        d.unc_theta = util::interp((*start).unc_theta, (*end).unc_theta, frac);
        d.rh = util::interp((*start).rh, (*end).rh, frac);
        d.delta = util::interp((*start).delta, (*end).delta, frac);
        d.phase = util::interp((*start).phase, (*end).phase, frac);
        d.selong = util::interp((*start).selong, (*end).selong, frac);
        d.true_anomaly = util::interp((*start).true_anomaly, (*end).true_anomaly, frac);
        d.sangle = util::interp((*start).sangle, (*end).sangle, frac);
        d.vangle = util::interp((*start).vangle, (*end).vangle, frac);
        d.vmag = util::interp((*start).vmag, (*end).vmag, frac);

        return Ephemeris{target_, {d}};
    }

    Ephemeris Ephemeris::extrapolate(const double distance, Ephemeris::Extrapolate direction) const
    {
        int i, j;
        if (direction == Ephemeris::Extrapolate::BACKWARDS)
        {
            i = 1;
            j = 0;
        }
        else
        {
            i = num_vertices() - 2;
            j = num_vertices() - 1;
        }
        const Datum d1 = data_[i];
        const Datum d2 = data_[j];
        const S2Point p1 = d1.as_s2point();
        const S2Point p2 = d2.as_s2point();
        const double length = S1Angle(p1, p2).radians();
        const double frac = 1 + distance / length;

        S2Point extrapolated = S2::Interpolate(p1, p2, frac).Normalize();

        Datum d;
        d.mjd = util::interp(d1.mjd, d2.mjd, frac);
        d.tmtp = util::interp(d1.tmtp, d2.tmtp, frac);
        d.radec(extrapolated);
        d.unc_a = util::interp(d1.unc_a, d2.unc_a, frac);
        d.unc_b = util::interp(d1.unc_b, d2.unc_b, frac);
        d.unc_theta = util::interp(d1.unc_theta, d2.unc_theta, frac);
        d.rh = util::interp(d1.rh, d2.rh, frac);
        d.delta = util::interp(d1.delta, d2.delta, frac);
        d.phase = util::interp(d1.phase, d2.phase, frac);
        d.selong = util::interp(d1.selong, d2.selong, frac);
        d.true_anomaly = util::interp(d1.true_anomaly, d2.true_anomaly, frac);
        d.sangle = util::interp(d1.sangle, d2.sangle, frac);
        d.vangle = util::interp(d1.vangle, d2.vangle, frac);
        d.vmag = util::interp(d1.vmag, d2.vmag, frac);

        return Ephemeris(target_, {d});
    }

    Ephemeris Ephemeris::subsample(const double mjd_start, const double mjd_stop) const
    {
        Ephemeris eph(target_, {});

        // find any whole segments between start and end
        vector<double> t = mjd();
        auto next = std::find_if(t.begin(), t.end(), [mjd_start](double t)
                                 { return (t >= mjd_start); });
        auto last = std::find_if(t.begin(), t.end(), [mjd_stop](double t)
                                 { return (t > mjd_stop); }) -
                    1;

        // Was interpolation between two epochs requested?
        if ((*next) > mjd_start)
            eph.append(interpolate(mjd_start));

        if (last >= next)
        {
            // there is at least one epoch between start and end
            for (int i = (next - t.begin()); i <= (last - t.begin()); i++)
                eph.append({data_[i]});
        }

        // Was interpolation between two epochs requested?
        if ((*last) < mjd_stop)
            eph.append(interpolate(mjd_stop));

        return eph;
    }

    void Ephemeris::pad(const vector<double> &a, const vector<double> &b, const vector<double> &theta, S2Polygon &polygon) const
    {
        // a, b in arcsec, theta in deg
        if ((a.size() != b.size()) | (a.size() != theta.size()) | (a.size() != num_vertices()))
            throw std::runtime_error("Length of padding vectors must match the length of the ephemeris (" +
                                     std::to_string(num_vertices_) +
                                     "): a, b, theta = " +
                                     std::to_string(a.size()) + " " +
                                     std::to_string(b.size()) + " " +
                                     std::to_string(theta.size()) + ".");

        const double max_a = *std::max_element(a.begin(), a.end());
        const double max_b = *std::max_element(b.begin(), b.end());
        if ((max_a * ARCSEC >= 90 * DEG) | (max_b * ARCSEC >= 90 * DEG))
            throw std::runtime_error("Padding must be less than 90 deg.");

        vector<S1Angle> theta_;
        for (double th : theta)
            theta_.push_back(S1Angle::Degrees(th));

        // Generate a convex hull around the error ellipse points
        S2ConvexHullQuery q;
        for (int i = 0; i < num_vertices_; i++)
        {
            for (auto p : util::ellipse(16, S2LatLng(data_[i].as_s2point()), a[i] * ARCSEC, b[i] * ARCSEC, theta[i] * DEG))
            {
                q.AddPoint(p.Normalized().ToPoint());
            }
        }

        auto loop = q.GetConvexHull();
        polygon.Init(std::move(loop));
    }

    void Ephemeris::pad(const double a, const double b, const double theta, S2Polygon &polygon) const
    {
        vector<double> a_vector(num_vertices_, a);
        vector<double> b_vector(num_vertices_, b);
        vector<double> theta_vector(num_vertices_, theta);
        pad(a_vector, b_vector, theta_vector, polygon);
    }

    void Ephemeris::pad(const vector<double> &para, const vector<double> &perp, S2Polygon &polygon) const
    {
        if ((para.size() != perp.size()) | (para.size() != num_vertices()))
            throw std::runtime_error("Length of padding vectors must match the length of the ephemeris.");

        vector<double> theta;
        for (int i = 0; i < para.size() - 1; i++)
            theta.push_back(util::position_angle(vertex(i), vertex(i + 1)) / DEG);
        // the PA of the last vertex is assumed to be the same as the one previous to it
        theta.push_back(*(theta.end() - 1));
        pad(para, perp, theta, polygon);
    }

    void Ephemeris::pad(const double para, const double perp, S2Polygon &polygon) const
    {
        vector<double> para_vector(num_vertices(), para);
        vector<double> perp_vector(num_vertices(), perp);
        pad(para_vector, perp_vector, polygon);
    }

    void Ephemeris::as_polygon(S2Polygon &polygon, vector<double> padding) const
    {
        // to generate a polygon, force the minimum padding to 0.1"
        vector<double> a(padding), b(padding), theta(num_vertices_, 0);

        std::replace_if(a.begin(), a.end(), [](double v)
                        { return v < 0.1; }, 0.1);
        std::replace_if(b.begin(), b.end(), [](double v)
                        { return v < 0.1; }, 0.1);

        if (options_.use_uncertainty)
        {
            // UNDEF_UNC evaluates to -1
            auto minimum_uncertainty = [this](const double x, const double y)
            { return ((x > y) ? x : y); };

            // make sure the values are >= 0.
            vector<double> unc = unc_a();
            std::transform(unc.begin(), unc.end(), a.begin(), a.begin(), minimum_uncertainty);

            unc = unc_b();
            std::transform(unc.begin(), unc.end(), b.begin(), b.begin(), minimum_uncertainty);
        }

        pad(a, b, theta, polygon);
    }

    void Ephemeris::as_polygon(S2Polygon &polygon) const
    {
        return as_polygon(polygon, vector<double>(num_vertices(), 0.1));
    }

    // void Ephemeris::as_polygons(vector<S2Polygon> &polygons) const
    // {
    // }

    json::array Ephemeris::as_json()
    {
        json::array data_array;
        for (Datum datum : data())
            data_array.emplace_back(datum.as_json());
        return data_array;
    }

    int Ephemeris::normalize_index(const int k, const int max) const
    {
        if ((k < -max) | (k >= max))
            throw std::runtime_error("Invalid index " + std::to_string(k) +
                                     " given number of elements: " + std::to_string(max));
        return k + ((k >= 0) ? 0 : max);
    }
}