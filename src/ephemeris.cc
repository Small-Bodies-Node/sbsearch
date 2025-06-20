#include "config.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <memory>
#include <optional>
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
#include <s2/s2region_union.h>
#include <s2/s2text_format.h>

#include "date.h"
#include "ephemeris.h"
#include "exceptions.h"
#include "observatory.h"
#include "table.h"
#include "util/math.h"
#include "util/spherical.h"
#include "util/string.h"

using sbsearch::table::Table;
using sbsearch::util::join;
using std::cerr;
using std::cout;
using std::endl;
using std::optional;
using std::unique_ptr;
using std::vector;

namespace sbsearch
{
    bool Ephemeris::Datum::operator==(const Ephemeris::Datum &other) const
    {
        return (std::tie(mjd, tmtp, ra, dec, mu, mu_theta,
                         unc_a, unc_b, unc_theta,
                         rh, delta, phase, selong, true_anomaly,
                         sangle, vangle, vmag) ==
                std::tie(other.mjd, other.tmtp, other.ra, other.dec, other.mu, other.mu_theta,
                         other.unc_a, other.unc_b, other.unc_theta,
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
        datum["mjd"] = json::value_from(mjd);
        datum["tmtp"] = json::value_from(tmtp);
        datum["ra"] = json::value_from(ra);
        datum["dec"] = json::value_from(dec);
        datum["mu"] = json::value_from(mu);
        datum["mu_theta"] = json::value_from(mu_theta);
        datum["unc_a"] = json::value_from(unc_a);
        datum["unc_b"] = json::value_from(unc_b);
        datum["unc_theta"] = json::value_from(unc_theta);
        datum["rh"] = json::value_from(rh);
        datum["delta"] = json::value_from(delta);
        datum["phase"] = json::value_from(phase);
        datum["selong"] = json::value_from(selong);
        datum["true_anomaly"] = json::value_from(true_anomaly);
        datum["sangle"] = json::value_from(sangle);
        datum["vangle"] = json::value_from(vangle);
        datum["vmag"] = json::value_from(vmag);
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

    std::istream &operator>>(std::istream &in, Ephemeris::Format::DateFormat &date_format)
    {
        std::string token;
        in >> token;
        std::transform(token.begin(), token.end(), token.begin(),
                       [](unsigned char c)
                       { return std::tolower(c); });
        if (token == "mjd")
            date_format = Ephemeris::Format::DateFormat::MJD;
        else if (token == "calendar")
            date_format = Ephemeris::Format::DateFormat::CALENDAR;
        else
            in.setstate(std::ios_base::failbit);
        return in;
    }

    std::ostream &operator<<(std::ostream &os, const Ephemeris &ephemeris)
    {
        Table table;
        if (ephemeris.format.date == Ephemeris::Format::DateFormat::CALENDAR)
            table.add_column("date", "%19s", ephemeris.date());
        else
            table.add_column("mjd", "%.6lf", ephemeris.mjd());
        table.add_column("tmtp", "%.6lf", ephemeris.tmtp());
        table.add_column("ra", "%.6lf", ephemeris.ra());
        table.add_column("dec", "%.6lf", ephemeris.dec());
        table.add_column("mu", "%.2f", ephemeris.mu());
        table.add_column("mu_theta", "%.3f", ephemeris.mu_theta());
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

    vector<optional<double>> Ephemeris::mjd() const
    {
        vector<optional<double>> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.mjd; });
        return result;
    }

    vector<string> Ephemeris::date() const
    {
        vector<string> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.mjd ? Date(datum.mjd.value()).iso() : "null"; });
        return result;
    }

    vector<optional<double>> Ephemeris::tmtp() const
    {
        vector<optional<double>> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.tmtp; });
        return result;
    }

    vector<optional<double>> Ephemeris::ra() const
    {
        vector<optional<double>> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.ra; });
        return result;
    }

    vector<optional<double>> Ephemeris::dec() const
    {
        vector<optional<double>> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.dec; });
        return result;
    }

    vector<optional<double>> Ephemeris::mu() const
    {
        vector<optional<double>> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.mu; });
        return result;
    }

    vector<optional<double>> Ephemeris::mu_theta() const
    {
        vector<optional<double>> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.mu_theta; });
        return result;
    }

    vector<optional<double>> Ephemeris::unc_a() const
    {
        vector<optional<double>> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.unc_a; });
        return result;
    }

    vector<optional<double>> Ephemeris::unc_b() const
    {
        vector<optional<double>> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.unc_b; });
        return result;
    }

    vector<optional<double>> Ephemeris::unc_theta() const
    {
        vector<optional<double>> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.unc_theta; });
        return result;
    }

    vector<optional<double>> Ephemeris::rh() const
    {
        vector<optional<double>> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.rh; });
        return result;
    }

    vector<optional<double>> Ephemeris::delta() const
    {
        vector<optional<double>> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.delta; });
        return result;
    }

    vector<optional<double>> Ephemeris::phase() const
    {
        vector<optional<double>> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.phase; });
        return result;
    }

    vector<optional<double>> Ephemeris::selong() const
    {
        vector<optional<double>> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.selong; });
        return result;
    }

    vector<optional<double>> Ephemeris::true_anomaly() const
    {
        vector<optional<double>> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.true_anomaly; });
        return result;
    }

    vector<optional<double>> Ephemeris::sangle() const
    {
        vector<optional<double>> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.sangle; });
        return result;
    }

    vector<optional<double>> Ephemeris::vangle() const
    {
        vector<optional<double>> result(num_vertices_);
        std::transform(data_.begin(), data_.end(), result.begin(),
                       [](Datum datum)
                       { return datum.vangle; });
        return result;
    }

    vector<optional<double>> Ephemeris::vmag() const
    {
        vector<optional<double>> result(num_vertices_);
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
        d.mu = util::interp((*start).mu, (*end).mu, frac);
        d.mu_theta = util::interp((*start).mu_theta, (*end).mu_theta, frac);
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
        d.mu = util::interp(d1.mu, d2.mu, frac);
        d.mu_theta = util::interp(d1.mu_theta, d2.mu_theta, frac);
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
        vector<optional<double>> t = mjd();
        auto next = std::lower_bound(t.begin(), t.end(), mjd_start);
        auto last = std::upper_bound(next, t.end(), mjd_stop) - 1;

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

    vector<unique_ptr<S2Polygon>> Ephemeris::as_polygons(double padding) const
    {
        // The minimum padding is a 0.1" radius circle
        padding = std::max(padding * 60.0, 0.1);
        vector<double> a(num_vertices_, padding), b(num_vertices_, padding), theta(num_vertices_, 0);

        if (options_.use_uncertainty)
        {
            for (int i = 0; i < num_vertices_; i++)
            {
                a[i] += data_[i].unc_a.value_or(0);
                b[i] += data_[i].unc_b.value_or(0);
                theta[i] = data_[i].unc_theta.value_or(0);
            }
        }

        // collect polygons in a vector
        std::vector<std::unique_ptr<S2Polygon>> polygons;

        // polygons around each ephemeris point
        for (int i = 0; i < num_vertices_; i++)
            polygons.emplace_back(
                util::circumscribe_ellipse(
                    S2LatLng(vertex(i)),
                    a[i] * ARCSEC,
                    b[i] * ARCSEC,
                    theta[i] * DEG,
                    data_[i].mu_theta.value() * DEG));

        // polygons between ephemeris points
        for (int i = 0; i < num_vertices_ - 1; i++)
        {
            std::vector<S2Point> vertices({{polygons[i]->loop(0)->vertex(0),
                                            polygons[i + 1]->loop(0)->vertex(1),
                                            polygons[i + 1]->loop(0)->vertex(2),
                                            polygons[i]->loop(0)->vertex(3)}});
            util::fix_crossing_edges(vertices);
            std::unique_ptr<S2Loop> loop = std::make_unique<S2Loop>(vertices, S2Debug::DISABLE);
            loop->Normalize();
            polygons.emplace_back(std::move(std::make_unique<S2Polygon>(S2Polygon(std::move(loop)))));
        }

        return std::move(polygons);
    }

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