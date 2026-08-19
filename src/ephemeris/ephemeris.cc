#include <iostream>
#include <optional>
#include <string>
#include <tuple>
#include <vector>
#include <boost/json.hpp>

#include "config.h"
#include "date.h"
#include "exceptions.h"
#include "moving_target.h"
#include "table.h"
#include "vertices.h"
#include "ephemeris/ephemeris.h"
#include "util/math.h"

using namespace sbsearch::table;

using std::optional;
using std::string;
using std::to_string;
using std::vector;

namespace sbsearch::ephemeris
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

    Ephemeris::Ephemeris(const MovingTarget &target, const Data &data, const Format f)
    {
        num_vertices_ = data.size();
        num_segments_ = (num_vertices_ == 0) ? 0 : (num_vertices_ - 1);

        target_ = target;
        data_ = data;
        format = f;

        isValid();
    }

    bool Ephemeris::isValid() const
    {
        if (!util::is_monotonically_increasing(mjd()))
            throw std::runtime_error("mjd must be monotonically increasing.");

        return true;
    }

    Ephemeris Ephemeris::operator[](const int k) const
    {
        return {target_, {data(k)}, format};
    }

    Ephemeris Ephemeris::slice(const int start) const
    {
        const int i = normalize_index(start, num_vertices_);
        Data subset(data_.begin() + i, data_.end());
        return {target_, subset, format};
    }

    Ephemeris Ephemeris::slice(const int start, const int stop) const
    {
        const int i = normalize_index(start, num_vertices_);
        const int j = normalize_index(stop, num_vertices_ + 1);

        if (i > j)
            throw EphemerisError("start cannot be greater than stop.");

        return {target_, {data_.begin() + i, data_.begin() + j}, format};
    }

    std::ostream &operator<<(std::ostream &os, const Ephemeris &ephemeris)
    {
        Table table;
        if (ephemeris.format.date == Date::Format::Calendar)
            table.add(Column("date", "%s", ephemeris.date()));
        else
            table.add(Column("mjd", "%.6lf", ephemeris.mjd()));
        table.add(Column("tmtp", "%.6lf", ephemeris.tmtp()));
        table.add(Column("ra", "%.6lf", ephemeris.ra()));
        table.add(Column("dec", "%.6lf", ephemeris.dec()));
        table.add(Column("mu", "%.2f", ephemeris.mu()));
        table.add(Column("mu_theta", "%.3f", ephemeris.mu_theta()));
        table.add(Column("rh", "%.4f", ephemeris.rh()));
        table.add(Column("delta", "%.4f", ephemeris.delta()));
        table.add(Column("phase", "%.3f", ephemeris.phase()));
        table.add(Column("selong", "%.3f", ephemeris.selong()));
        table.add(Column("true_anomaly", "%.3f", ephemeris.true_anomaly()));
        table.add(Column("sangle", "%.3f", ephemeris.sangle()));
        table.add(Column("vangle", "%.3f", ephemeris.vangle()));
        table.add(Column("unc_a", "%.3f", ephemeris.unc_a()));
        table.add(Column("unc_b", "%.3f", ephemeris.unc_b()));
        table.add(Column("unc_th", "%.3f", ephemeris.unc_theta()));
        table.add(Column("vmag", "%.3f", ephemeris.vmag()));

        os << table;
        return os;
    }

    bool Ephemeris::operator==(const Ephemeris &other) const
    {
        return ((target_ == other.target()) && (data_ == other.data()));
    }

    const Ephemeris::Datum &Ephemeris::data(const int k) const
    {
        const int i = normalize_index(k, num_vertices_);
        return data_[i];
    }

    vector<optional<double>> Ephemeris::mjd() const
    {
        return get_data_vector<vector<optional<double>>>(&Datum::mjd);
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
        return get_data_vector<vector<optional<double>>>(&Datum::tmtp);
    }

    vector<optional<double>> Ephemeris::ra() const
    {
        return get_data_vector<vector<optional<double>>>(&Datum::ra);
    }

    vector<optional<double>> Ephemeris::dec() const
    {
        return get_data_vector<vector<optional<double>>>(&Datum::dec);
    }

    vector<optional<double>> Ephemeris::mu() const
    {
        return get_data_vector<vector<optional<double>>>(&Datum::mu);
    }

    vector<optional<double>> Ephemeris::mu_theta() const
    {
        return get_data_vector<vector<optional<double>>>(&Datum::mu_theta);
    }

    vector<optional<double>> Ephemeris::unc_a() const
    {
        return get_data_vector<vector<optional<double>>>(&Datum::unc_a);
    }

    vector<optional<double>> Ephemeris::unc_b() const
    {
        return get_data_vector<vector<optional<double>>>(&Datum::unc_b);
    }

    vector<optional<double>> Ephemeris::unc_theta() const
    {
        return get_data_vector<vector<optional<double>>>(&Datum::unc_theta);
    }

    vector<optional<double>> Ephemeris::rh() const
    {
        return get_data_vector<vector<optional<double>>>(&Datum::rh);
    }

    vector<optional<double>> Ephemeris::delta() const
    {
        return get_data_vector<vector<optional<double>>>(&Datum::delta);
    }

    vector<optional<double>> Ephemeris::phase() const
    {
        return get_data_vector<vector<optional<double>>>(&Datum::phase);
    }

    vector<optional<double>> Ephemeris::selong() const
    {
        return get_data_vector<vector<optional<double>>>(&Datum::selong);
    }

    vector<optional<double>> Ephemeris::true_anomaly() const
    {
        return get_data_vector<vector<optional<double>>>(&Datum::true_anomaly);
    }

    vector<optional<double>> Ephemeris::sangle() const
    {
        return get_data_vector<vector<optional<double>>>(&Datum::sangle);
    }

    vector<optional<double>> Ephemeris::vangle() const
    {
        return get_data_vector<vector<optional<double>>>(&Datum::vangle);
    }

    vector<optional<double>> Ephemeris::vmag() const
    {
        return get_data_vector<vector<optional<double>>>(&Datum::vmag);
    }

    int Ephemeris::num_vertices() const
    {
        return num_vertices_;
    }

    S2Point Ephemeris::vertex(const int k) const
    {
        return make_point(data(k));
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

    Ephemeris Ephemeris::segment(const int k) const
    {
        const int i = normalize_index(k, num_segments_);
        return {target_, {data_[i], data_[i + 1]}, format};
    }

    vector<Ephemeris> Ephemeris::segments() const
    {
        vector<Ephemeris> eph(num_segments_);
        for (int i = 0; i < num_segments_; i++)
            eph[i] = segment(i);
        return eph;
    }

    void Ephemeris::append(const Datum &new_datum)
    {
        check_relative_time(*this, {new_datum});
        data_.push_back(new_datum);
        num_vertices_ = data_.size();
        num_segments_ = (num_vertices_ == 0) ? 0 : (num_vertices_ - 1);
    }

    void Ephemeris::append(Datum &&new_datum)
    {
        check_relative_time(*this, {new_datum});
        data_.push_back(std::move(new_datum));
        num_vertices_ = data_.size();
        num_segments_ = (num_vertices_ == 0) ? 0 : (num_vertices_ - 1);
    }

    void Ephemeris::append(const Data &new_data)
    {
        check_relative_time(*this, new_data);
        std::copy(new_data.begin(), new_data.end(), std::back_inserter(data_));
        num_vertices_ = data_.size();
        num_segments_ = (num_vertices_ == 0) ? 0 : (num_vertices_ - 1);
    }

    void Ephemeris::append(Data &&new_data)
    {
        check_relative_time(*this, new_data);
        std::move(new_data.begin(), new_data.end(), std::back_inserter(data_));
        num_vertices_ = data_.size();
        num_segments_ = (num_vertices_ == 0) ? 0 : (num_vertices_ - 1);
    }

    void Ephemeris::append(const Ephemeris &eph)
    {
        check_target_id(*this, eph);
        append(eph.data());
    }

    void Ephemeris::append(Ephemeris &&eph)
    {
        check_target_id(*this, eph);
        append(std::move(eph.data()));
    }

    int Ephemeris::normalize_index(const int k, const int max)
    {
        if ((k < -max) || (k >= max))
            throw std::runtime_error("Invalid index " + std::to_string(k) +
                                     " given number of elements: " + std::to_string(max));
        return k + ((k >= 0) ? 0 : max);
    }

    void Ephemeris::check_target_id(const Ephemeris &a, const Ephemeris &b)
    {
        if (a.target().moving_target_id() != b.target().moving_target_id())
            throw std::runtime_error("Attempted to append an ephemeris with a different object ID: " +
                                     a.target().to_string() + " and " + b.target().to_string());
    }

    void Ephemeris::check_relative_time(const Ephemeris &a, const Ephemeris::Data &b)
    {
        // check that b's time axis follows a
        if (a.num_vertices() > 0 && (a.data().back().mjd > b[0].mjd))
            throw std::runtime_error("Attempting to append ephemeris data with an earlier mjd. Compare" +
                                     to_string(a.data().back().mjd.value()) + " to " +
                                     to_string(b[0].mjd.value()));

        // check that b's time axis is in order
        if (b.size() > 1)
        {
            auto i = std::adjacent_find(b.begin(), b.end(),
                                        [](const auto &left, const auto &right)
                                        { return left.mjd > right.mjd; });
            if (i != b.end())
                throw std::runtime_error("mjd must be monotonically increasing.");
        }
    }

}