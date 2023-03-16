#include "config.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <stdexcept>
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
#include "sbsearch_testing.h"
#include "util.h"

using sbsearch::position_angle;
using std::cerr;
using std::cout;
using std::endl;
using std::make_unique;
using std::unique_ptr;
using std::vector;

namespace sbsearch
{
    Ephemeris::Ephemeris(const int object_id,
                         const vector<S2Point> &vertices,
                         const vector<double> &mjd,
                         const vector<double> &rh,
                         const vector<double> &delta,
                         const vector<double> &phase,
                         const vector<double> &unc_a,
                         const vector<double> &unc_b,
                         const vector<double> &unc_theta)
    {
        num_vertices_ = vertices.size();
        num_segments_ = (num_vertices_ == 0) ? 0 : (num_vertices_ - 1);

        if ((mjd.size() != num_vertices_) |
            (rh.size() != num_vertices_) |
            (delta.size() != num_vertices_) |
            (phase.size() != num_vertices_) |
            (unc_a.size() != num_vertices_) |
            (unc_b.size() != num_vertices_) |
            (unc_theta.size() != num_vertices_))
            throw std::runtime_error("The size of all vectors must match.");

        object_id_ = object_id;
        vertices_ = vector<S2Point>(vertices.begin(), vertices.end());
        mjd_ = vector<double>(mjd.begin(), mjd.end());
        rh_ = vector<double>(rh.begin(), rh.end());
        delta_ = vector<double>(delta.begin(), delta.end());
        phase_ = vector<double>(phase.begin(), phase.end());
        unc_a_ = vector<double>(unc_a.begin(), unc_a.end());
        unc_b_ = vector<double>(unc_b.begin(), unc_b.end());
        unc_theta_ = vector<double>(unc_theta.begin(), unc_theta.end());

        isValid();
    }

    Ephemeris Ephemeris::operator[](const int k) const
    {
        if ((k < -num_vertices()) | (k >= num_vertices()))
            throw std::runtime_error("Invalid index.");

        const int i = k + ((k >= 0) ? 0 : num_vertices());
        Ephemeris eph = Ephemeris(object_id_,
                                  {vertices_[i]},
                                  {mjd_[i]},
                                  {rh_[i]},
                                  {delta_[i]},
                                  {phase_[i]},
                                  {unc_a_[i]},
                                  {unc_b_[i]},
                                  {unc_theta_[i]});
        eph.format = format;
        return eph;
    }

    bool Ephemeris::isValid() const
    {
        if (!sbsearch::is_increasing(mjd_))
            throw std::runtime_error("mjd must be monotonically increasing.");

        return true;
    }

    std::ostream &operator<<(std::ostream &os, const Ephemeris &ephemeris)
    {
        if (ephemeris.num_vertices() == 1)
        {
            os << std::fixed
               << std::right
               << std::setw(ephemeris.format.object_id_width)
               << ephemeris.object_id() << "  "
               << std::setw(11)
               << std::setprecision(5)
               << ephemeris.mjd(0) << "  "
               << std::setw(12)
               << std::setprecision(6)
               << ephemeris.ra(0) << "  "
               << std::setw(12)
               << ephemeris.dec(0) << "  "
               << std::setw(6)
               << std::setprecision(3)
               << ephemeris.rh(0) << "  "
               << std::setw(6)
               << ephemeris.delta(0) << "  "
               << std::setw(6)
               << std::setprecision(2)
               << ephemeris.phase(0)
               << std::defaultfloat;
        }
        else
        {
            for (int i = 0; i < ephemeris.num_vertices(); i++)
                os << ephemeris[i] << "\n";
        }

        return os;
    }
    bool Ephemeris::is_equal(const Ephemeris &other) const
    {
        return ((object_id() == other.object_id()) &
                (num_vertices() == other.num_vertices()) &
                (vertices() == other.vertices()) &
                (mjd() == other.mjd()) &
                (rh() == other.rh()) &
                (delta() == other.delta()) &
                (phase() == other.phase()) &
                (unc_a() == other.unc_a()) &
                (unc_b() == other.unc_b()) &
                (unc_theta() == other.unc_theta()));
    }

    int Ephemeris::num_vertices() const
    {
        return num_vertices_;
    }

    const S2Point &Ephemeris::vertex(const int k) const
    {
        if ((k < -num_vertices()) | (k >= num_vertices()))
            throw std::runtime_error("Invalid index.");

        return vertices_[k + ((k >= 0) ? 0 : num_vertices())];
    }

    double Ephemeris::ra(const int k) const
    {
        // S2 is -180 to 180, but we like 0 to 360
        return fmod(S2LatLng::Longitude(vertex(k)).degrees() + 360, 360);
    }

    double Ephemeris::dec(const int k) const
    {
        return S2LatLng::Latitude(vertex(k)).degrees();
    }

    int Ephemeris::num_segments() const
    {
        return num_segments_;
    }

    void Ephemeris::append(const Ephemeris &eph)
    {
        if (eph.object_id() != object_id_)
            throw std::runtime_error("Attempted to append an ephemeris with a different object ID.");

        if (num_vertices_ != 0)
        {
            if (mjd(-1) > eph.mjd(0))
                throw std::runtime_error("Attempting to append an ephemeris with an earlier mjd.");
        }

        vertices_.insert(vertices_.end(), eph.vertices().begin(), eph.vertices().end());
        mjd_.insert(mjd_.end(), eph.mjd().begin(), eph.mjd().end());
        rh_.insert(rh_.end(), eph.rh().begin(), eph.rh().end());
        delta_.insert(delta_.end(), eph.delta().begin(), eph.delta().end());
        phase_.insert(phase_.end(), eph.phase().begin(), eph.phase().end());
        unc_a_.insert(unc_a_.end(), eph.unc_a().begin(), eph.unc_a().end());
        unc_b_.insert(unc_b_.end(), eph.unc_b().begin(), eph.unc_b().end());
        unc_theta_.insert(unc_theta_.end(), eph.unc_theta().begin(), eph.unc_theta().end());

        num_vertices_ = vertices_.size();
        num_segments_ = (num_vertices_ == 0) ? 0 : (num_vertices_ - 1);
    }

    Ephemeris Ephemeris::segment(const int k) const
    {
        if ((k < -num_segments()) | (k >= num_segments()))
            throw std::runtime_error("Invalid index.");

        const int i = k + ((k < 0) ? num_segments() : 0);
        Ephemeris eph = Ephemeris(object_id_,
                                  {vertices_[i], vertices_[i + 1]},
                                  {mjd_[i], mjd_[i + 1]},
                                  {rh_[i], rh_[i + 1]},
                                  {delta_[i], delta_[i + 1]},
                                  {phase_[i], phase_[i + 1]},
                                  {unc_a_[i], unc_a_[i + 1]},
                                  {unc_b_[i], unc_b_[i + 1]},
                                  {unc_theta_[i], unc_theta_[i + 1]});
        eph.format = format;
        return eph;
    }

    vector<Ephemeris> Ephemeris::segments() const
    {
        vector<Ephemeris> eph;
        for (int i = 0; i < num_segments(); i++)
            eph.push_back(segment(i));
        return eph;
    }

    S2Polyline Ephemeris::as_polyline() const
    {
        return S2Polyline(vertices());
    }

    Ephemeris Ephemeris::interpolate(const double mjd0) const
    {
        if ((mjd0 < *mjd_.begin()) | (mjd0 > *(mjd_.end() - 1)))
            throw std::runtime_error("Interpolation beyond ephemeris time range.");

        // find the nearest segment
        int start = std::find_if(mjd_.begin(), mjd_.end(), [mjd0](double t)
                                 { return (t > mjd0); }) -
                    1 - mjd_.begin();
        int end = start + 1;

        // length of segment in time
        double dt = mjd_[end] - mjd_[start];

        // interpolate to this fraction
        double frac = (mjd0 - mjd_[start]) / dt;

        // this is the line we will interpolate
        S2Polyline segment(vector<S2Point>(vertices_.begin() + start, vertices_.begin() + end + 1));
        Ephemeris eph = Ephemeris(object_id_,
                                  {segment.Interpolate(frac)},
                                  {interp(mjd_[start], mjd_[end], frac)},
                                  {interp(rh_[start], rh_[end], frac)},
                                  {interp(delta_[start], delta_[end], frac)},
                                  {interp(phase_[start], phase_[end], frac)},
                                  {interp(unc_a_[start], unc_a_[end], frac)},
                                  {interp(unc_b_[start], unc_b_[end], frac)},
                                  {interp(unc_theta_[start], unc_theta_[end], frac)});
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
        const double length = S1Angle(vertices_[i], vertices_[j]).radians();
        const double frac = 1 + distance / length;

        S2Point extrapolated = S2::Interpolate(vertices_[i], vertices_[j], frac).Normalize();

        Ephemeris eph = Ephemeris(object_id_,
                                  {extrapolated},
                                  {interp(mjd_[i], mjd_[j], frac)},
                                  {interp(rh_[i], rh_[j], frac)},
                                  {interp(delta_[i], delta_[j], frac)},
                                  {interp(phase_[i], phase_[j], frac)},
                                  {interp(unc_a_[i], unc_a_[j], frac)},
                                  {interp(unc_b_[i], unc_b_[j], frac)},
                                  {interp(unc_theta_[i], unc_theta_[j], frac)});
        eph.format = format;
        return eph;
    }

    Ephemeris Ephemeris::subsample(const double mjd_start, const double mjd_stop) const
    {
        Ephemeris eph;
        eph.object_id(object_id_);
        eph.format = format;

        // find any whole segments between start and end
        auto next = std::find_if(mjd_.begin(), mjd_.end(), [mjd_start](double t)
                                 { return (t >= mjd_start); });
        auto last = std::find_if(mjd_.begin(), mjd_.end(), [mjd_stop](double t)
                                 { return (t > mjd_stop); }) -
                    1;

        // Was interpolation between two epochs requested?
        if ((*next) > mjd_start)
            eph.append(interpolate(mjd_start));

        if (last >= next)
        {
            // there is at least one epoch between start and end
            for (int i = (next - mjd_.begin()); i <= (last - mjd_.begin()); i++)
                eph.append((*this)[i]);
        }

        // Was interpolation between two epochs requested?
        if ((*last) < mjd_stop)
            eph.append(interpolate(mjd_stop));

        return eph;
    }

    S2Polygon Ephemeris::pad(const vector<double> &a, const vector<double> &b, const vector<double> &theta) const
    {
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
            for (auto p : ellipse(16, S2LatLng(vertices_[i]), a[i] * ARCSEC, b[i] * ARCSEC, theta[i] * DEG))
            {
                q.AddPoint(p.Normalized().ToPoint());
            }
        }

        auto loop = q.GetConvexHull();
        S2Polygon polygon(std::move(loop));
        return polygon;
    }

    S2Polygon Ephemeris::pad(const double a, const double b, const double theta) const
    {
        vector<double> a_vector(num_vertices(), a);
        vector<double> b_vector(num_vertices(), b);
        vector<double> theta_vector(num_vertices(), theta);
        return pad(a_vector, b_vector, theta_vector);
    }

    S2Polygon Ephemeris::pad(const vector<double> &para, const vector<double> &perp) const
    {
        if ((para.size() != perp.size()) | (para.size() != num_vertices()))
            throw std::runtime_error("Length of padding vectors must match the length of the ephemeris.");

        vector<double> theta;
        for (int i = 0; i < para.size() - 1; i++)
            theta.push_back(position_angle(vertex(i), vertex(i + 1)) / DEG);
        // the PA of the last vertex is assumed to be the same as the one previous to it
        theta.push_back(*(theta.end() - 1));
        return pad(para, perp, theta);
    }

    S2Polygon Ephemeris::pad(const double para, const double perp) const
    {
        vector<double> para_vector(num_vertices(), para);
        vector<double> perp_vector(num_vertices(), perp);
        return pad(para_vector, perp_vector);
    }

    S2Polygon Ephemeris::as_polygon() const
    {
        // minimum padding is 0.1 arcsec

        S2Polygon *polygon = new S2Polygon;
        vector<double> a(num_vertices_, 0.1), b(num_vertices_, 0.1), theta(num_vertices_, 0);

        if (options_.use_uncertainty)
        {
            auto add_uncertainty = [this](const double x, const double y)
            { return ((x >= 0) ? x : 0) + y; };

            // make sure the values are >= 0.
            std::transform(unc_a_.begin(), unc_a_.end(), a.begin(), a.begin(), add_uncertainty);
            std::transform(unc_b_.begin(), unc_b_.end(), b.begin(), b.begin(), add_uncertainty);
            std::copy(unc_theta_.begin(), unc_theta_.end(), theta.begin());
        }

        if (options_.padding > 0)
        {
            auto add_padding = [this](const double x)
            { return x + this->options_.padding; };

            std::transform(a.begin(), a.end(), a.begin(), add_padding);
            std::transform(b.begin(), b.end(), b.begin(), add_padding);
        }

        return pad(a, b, theta);
    }

    const double &Ephemeris::getter(const vector<double> &property, const int k) const
    {
        if ((k < -num_vertices()) | (k >= num_vertices()))
            throw std::runtime_error("Invalid index.");
        return property[k + ((k >= 0) ? 0 : num_vertices())];
    }
}