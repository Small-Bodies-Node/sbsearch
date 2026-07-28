#include "config.h"
#include <optional>
#include "ephemeris.h"

using std::optional;
using std::vector;

namespace sbsearch
{
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

}