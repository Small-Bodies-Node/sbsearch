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
#include "util.h"

using sbsearch::position_angle;
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

    const Ephemeris Ephemeris::slice(const int start)
    {
        const int i = normalize_index(start, num_vertices_);
        Data subset(data_.begin() + i, data_.end());
        return Ephemeris(target_, subset);
    }

    const Ephemeris Ephemeris::slice(const int start, const int stop)
    {
        const int i = normalize_index(start, num_vertices_);
        const int j = normalize_index(stop, num_vertices_);

        if (i > j)
            throw EphemerisError("start cannot be greater than stop.");

        Data subset(data_.begin() + i, data_.begin() + j);
        return Ephemeris(target_, subset);
    }

    bool Ephemeris::isValid() const
    {
        if (!sbsearch::is_increasing(mjd()))
            throw std::runtime_error("mjd must be monotonically increasing.");

        return true;
    }

    Ephemeris::Format Ephemeris::format_widths() const
    {
        vector<double> abs_tmtp = tmtp();
        std::transform(abs_tmtp.begin(), abs_tmtp.end(), abs_tmtp.begin(), [](double x)
                       { return std::fabs(x); });

        Format format_{
            target_.designation().size(),
            std::to_string(target_.moving_target_id()).size(),
            (size_t)(std::floor(std::log10(*std::max_element(abs_tmtp.begin(), abs_tmtp.end()))) + 8),
        };

        return format_;
    }

    std::ostream &operator<<(std::ostream &os, const Ephemeris &ephemeris)
    {
        for (int i = 0; i < ephemeris.num_vertices(); i++)
        {
            Ephemeris::Datum row = ephemeris.data(i);
            os << std::fixed
               << std::right
               << std::setw(ephemeris.format.designation_width)
               << ephemeris.target().designation() << "  "
               << std::setw(ephemeris.format.moving_target_id_width)
               << ephemeris.target().moving_target_id() << "  "
               << std::setw(11)
               << std::setprecision(5)
               << row.mjd << "  "
               << std::setw(ephemeris.format.tmtp_width)
               << std::setprecision(5)
               << row.tmtp << "  "
               << std::setw(10)
               << std::setprecision(6)
               << row.ra << "  "
               << std::setw(10)
               << row.dec << "  "
               << std::setw(6)
               << std::setprecision(3)
               << row.rh << "  "
               << std::setw(6)
               << row.delta << "  "
               << std::setw(6)
               << std::setprecision(2)
               << row.phase
               << std::defaultfloat;

            if (ephemeris.num_vertices() != 1)
                os << "\n";
        }

        return os;
    }

    bool Ephemeris::operator==(const Ephemeris &other) const
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
                       { return datum.mjd; });
        return result;
    }

    vector<double> Ephemeris::tmtp() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.tmtp; });
        return result;
    }

    vector<double> Ephemeris::ra() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.ra; });
        return result;
    }

    vector<double> Ephemeris::dec() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.dec; });
        return result;
    }

    vector<double> Ephemeris::unc_a() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.unc_a; });
        return result;
    }

    vector<double> Ephemeris::unc_b() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.unc_b; });
        return result;
    }

    vector<double> Ephemeris::unc_theta() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.unc_theta; });
        return result;
    }

    vector<double> Ephemeris::rh() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.rh; });
        return result;
    }

    vector<double> Ephemeris::delta() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.delta; });
        return result;
    }

    vector<double> Ephemeris::phase() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.phase; });
        return result;
    }

    vector<double> Ephemeris::selong() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.selong; });
        return result;
    }

    vector<double> Ephemeris::true_anomaly() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.true_anomaly; });
        return result;
    }

    vector<double> Ephemeris::sangle() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.sangle; });
        return result;
    }

    vector<double> Ephemeris::vangle() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.vangle; });
        return result;
    }

    vector<double> Ephemeris::vmag() const
    {
        vector<double> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.vmag; });
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
        Ephemeris eph = Ephemeris(target_, {data_[i], data_[i + 1]});
        eph.format = format;
        return eph;
    }

    vector<Ephemeris> Ephemeris::segments() const
    {
        vector<Ephemeris> eph(num_segments_);
        for (int i = 0; i < num_segments_; i++)
            eph[i] = segment(i);
        return eph;
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
                new_data[i].mjd,
                new_data[i].delta);
            new_data[i].radec(coords.Normalized().ToPoint());
        }
        return Ephemeris(target_, new_data);
    }

    Ephemeris Ephemeris::interpolate(const double mjd) const
    {
        if ((mjd < data_.front().mjd) | (mjd > data_.back().mjd))
            throw std::runtime_error("Interpolation beyond ephemeris time range.");

        // find the nearest segment
        auto end = std::find_if(data_.begin(), data_.end(), [mjd](Datum d)
                                { return (d.mjd > mjd); });
        auto start = end - 1;

        // length of segment in time
        double dt = (*end).mjd - (*start).mjd;

        // interpolate to this fraction
        double frac = (mjd - (*start).mjd) / dt;

        // this is the line we will interpolate
        S2Polyline segment(vector<S2Point>({(*start).as_s2point(), (*end).as_s2point()}));

        Datum d;
        d.mjd = interp((*start).mjd, (*end).mjd, frac);
        d.tmtp = interp((*start).tmtp, (*end).tmtp, frac);
        d.radec(segment.Interpolate(frac));
        d.unc_a = interp((*start).unc_a, (*end).unc_a, frac);
        d.unc_b = interp((*start).unc_b, (*end).unc_b, frac);
        d.unc_theta = interp((*start).unc_theta, (*end).unc_theta, frac);
        d.rh = interp((*start).rh, (*end).rh, frac);
        d.delta = interp((*start).delta, (*end).delta, frac);
        d.phase = interp((*start).phase, (*end).phase, frac);
        d.selong = interp((*start).selong, (*end).selong, frac);
        d.true_anomaly = interp((*start).true_anomaly, (*end).true_anomaly, frac);
        d.sangle = interp((*start).sangle, (*end).sangle, frac);
        d.vangle = interp((*start).vangle, (*end).vangle, frac);
        d.vmag = interp((*start).vmag, (*end).vmag, frac);

        Ephemeris eph{target_, {d}};
        eph.format = format;
        return eph;
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
        d.mjd = interp(d1.mjd, d2.mjd, frac);
        d.tmtp = interp(d1.tmtp, d2.tmtp, frac);
        d.radec(extrapolated);
        d.unc_a = interp(d1.unc_a, d2.unc_a, frac);
        d.unc_b = interp(d1.unc_b, d2.unc_b, frac);
        d.unc_theta = interp(d1.unc_theta, d2.unc_theta, frac);
        d.rh = interp(d1.rh, d2.rh, frac);
        d.delta = interp(d1.delta, d2.delta, frac);
        d.phase = interp(d1.phase, d2.phase, frac);
        d.selong = interp(d1.selong, d2.selong, frac);
        d.true_anomaly = interp(d1.true_anomaly, d2.true_anomaly, frac);
        d.sangle = interp(d1.sangle, d2.sangle, frac);
        d.vangle = interp(d1.vangle, d2.vangle, frac);
        d.vmag = interp(d1.vmag, d2.vmag, frac);

        Ephemeris eph = Ephemeris(target_, {d});
        eph.format = format;
        return eph;
    }

    Ephemeris Ephemeris::subsample(const double mjd_start, const double mjd_stop) const
    {
        Ephemeris eph(target_, {});
        eph.format = format;

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
            throw std::runtime_error("Length of padding vectors must match the length of the ephemeris.");

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
            for (auto p : ellipse(16, S2LatLng(data_[i].as_s2point()), a[i] * ARCSEC, b[i] * ARCSEC, theta[i] * DEG))
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
            theta.push_back(position_angle(vertex(i), vertex(i + 1)) / DEG);
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

    void Ephemeris::as_polygon(S2Polygon &polygon) const
    {
        // to generate a polygon, force the minimum padding to 0.1"
        vector<double> a(num_vertices_, 0.1);
        vector<double> b(num_vertices_, 0.1);
        vector<double> theta(num_vertices_, 0);

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

        if (options_.padding > 0.1)
        {
            auto add_padding = [this](const double x)
            { return x + this->options_.padding; };

            std::transform(a.begin(), a.end(), a.begin(), add_padding);
            std::transform(b.begin(), b.end(), b.begin(), add_padding);
        }

        pad(a, b, theta, polygon);
    }

    int Ephemeris::normalize_index(const int k, const int max) const
    {
        if ((k < -max) | (k >= max))
            throw std::runtime_error("Invalid index.");
        return k + ((k >= 0) ? 0 : max);
    }

}