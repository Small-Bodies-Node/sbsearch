#include "config.h"

#include <cinttypes>
#include <optional>
#include <ostream>
#include <boost/json.hpp>

#include "ephemeris.h"
#include "found.h"
#include "observation.h"
#include "table.h"

using sbsearch::table::Table;
using std::optional;

namespace sbsearch
{
    const bool Found::operator==(const Found &other) const
    {
        return ((observation == other.observation) & (ephemeris == other.ephemeris));
    };

    const bool Found::operator!=(const Found &other) const
    {
        return !(*this == other);
    };

    json::object Found::as_json()
    {
        json::object obj;

        // if found.ephemeris is a segment, interpolate it to observation mid-time.
        Ephemeris eph;
        if (ephemeris.num_vertices() > 1)
        {
            double mjd = (observation.mjd_start() + observation.mjd_stop()) / 2;
            eph = ephemeris.interpolate(mjd);
        }
        else
            eph = ephemeris;

        for (auto item : observation.as_json())
            obj[item.key()] = item.value();

        json::object eph_object = *eph.as_json().at(0).if_object();
        for (auto item : eph_object)
            obj[item.key()] = item.value();

        return obj;
    }

    Founds::Founds(const vector<Found> &founds)
    {
        append(founds);
    }

    Founds::Founds(const Founds &founds)
    {
        append(founds.data);
    }

    void Founds::append(const Found &found)
    {
        data.push_back(found);
    }

    void Founds::append(const vector<Found> &founds)
    {
        data.reserve(data.size() + founds.size());
        data.insert(data.end(), founds.begin(), founds.end());
    }

    void Founds::append(const Founds &founds)
    {
        append(founds.data);
    }

    vector<string> Founds::source() const
    {
        int n = data.size();
        vector<string> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       { return found.observation.source(); });
        return v;
    }

    vector<string> Founds::observatory() const
    {
        int n = data.size();
        vector<string> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       { return found.observation.observatory(); });
        return v;
    }

    vector<string> Founds::product_id() const
    {
        int n = data.size();
        vector<string> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       { return found.observation.product_id(); });
        return v;
    }

    vector<string> Founds::fov() const
    {
        int n = data.size();
        vector<string> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       { return found.observation.fov(); });
        return v;
    }

    vector<optional<int64_t>> Founds::observation_id() const
    {
        int n = data.size();
        vector<optional<int64_t>> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       { return found.observation.observation_id(); });

        return v;
    }

    vector<double> Founds::mjd_start() const
    {
        int n = data.size();
        vector<double> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       { return found.observation.mjd_start(); });
        return v;
    }

    vector<double> Founds::mjd_stop() const
    {
        int n = data.size();
        vector<double> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       { return found.observation.mjd_stop(); });
        return v;
    }

    vector<double> Founds::exposure() const
    {
        int n = data.size();
        vector<double> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       { return found.observation.exposure(); });

        return v;
    }

    vector<optional<int64_t>> Founds::moving_target_id() const
    {
        int n = data.size();
        vector<optional<int64_t>> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       { return found.ephemeris.target().moving_target_id(); });
        return v;
    }

    vector<string> Founds::designation() const
    {
        int n = data.size();
        vector<string> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       { return found.ephemeris.target().designation(); });
        return v;
    }

    vector<bool> Founds::small_body() const
    {
        int n = data.size();
        vector<bool> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       { return found.ephemeris.target().small_body(); });
        return v;
    }

    vector<optional<double>> Founds::mjd() const
    {
        int n = data.size();
        vector<optional<double>> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       {
                           Ephemeris eph = (found.ephemeris.num_vertices() == 1)
                                               ? found.ephemeris
                                               : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return eph.data(0).mjd.value_or(-1);
                       });
        return v;
    }

    vector<optional<double>> Founds::tmtp() const
    {
        int n = data.size();
        vector<optional<double>> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       {
                           Ephemeris eph = (found.ephemeris.num_vertices() == 1)
                                               ? found.ephemeris
                                               : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return eph.data(0).tmtp.value_or(-1);
                       });
        return v;
    }

    vector<optional<double>> Founds::ra() const
    {
        int n = data.size();
        vector<optional<double>> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       {
                           Ephemeris eph = (found.ephemeris.num_vertices() == 1)
                                               ? found.ephemeris
                                               : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return eph.data(0).ra;
                       });
        return v;
    }

    vector<optional<double>> Founds::dec() const
    {
        int n = data.size();
        vector<optional<double>> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       {
                           Ephemeris eph = (found.ephemeris.num_vertices() == 1)
                                               ? found.ephemeris
                                               : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return eph.data(0).dec;
                       });
        return v;
    }

    vector<optional<double>> Founds::mu() const
    {
        int n = data.size();
        vector<optional<double>> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       {
                           Ephemeris eph = (found.ephemeris.num_vertices() == 1)
                                               ? found.ephemeris
                                               : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return eph.data(0).mu;
                       });
        return v;
    }

    vector<optional<double>> Founds::mu_theta() const
    {
        int n = data.size();
        vector<optional<double>> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       {
                           Ephemeris eph = (found.ephemeris.num_vertices() == 1)
                                               ? found.ephemeris
                                               : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return eph.data(0).mu_theta;
                       });
        return v;
    }

    vector<optional<double>> Founds::unc_a() const
    {
        int n = data.size();
        vector<optional<double>> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       {
                           Ephemeris eph = (found.ephemeris.num_vertices() == 1)
                                               ? found.ephemeris
                                               : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return eph.data(0).unc_a;
                       });
        return v;
    }

    vector<optional<double>> Founds::unc_b() const
    {
        int n = data.size();
        vector<optional<double>> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       {
                           Ephemeris eph = (found.ephemeris.num_vertices() == 1)
                                               ? found.ephemeris
                                               : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return eph.data(0).unc_b;
                       });
        return v;
    }

    vector<optional<double>> Founds::unc_theta() const
    {
        int n = data.size();
        vector<optional<double>> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       {
                           Ephemeris eph = (found.ephemeris.num_vertices() == 1)
                                               ? found.ephemeris
                                               : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return eph.data(0).unc_theta;
                       });
        return v;
    }

    vector<optional<double>> Founds::rh() const
    {
        int n = data.size();
        vector<optional<double>> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       {
                           Ephemeris eph = (found.ephemeris.num_vertices() == 1)
                                               ? found.ephemeris
                                               : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return eph.data(0).rh;
                       });
        return v;
    }

    vector<optional<double>> Founds::delta() const
    {
        int n = data.size();
        vector<optional<double>> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       {
                           Ephemeris eph = (found.ephemeris.num_vertices() == 1)
                                               ? found.ephemeris
                                               : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return eph.data(0).delta;
                       });
        return v;
    }

    vector<optional<double>> Founds::phase() const
    {
        int n = data.size();
        vector<optional<double>> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       {
                           Ephemeris eph = (found.ephemeris.num_vertices() == 1)
                                               ? found.ephemeris
                                               : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return eph.data(0).phase;
                       });
        return v;
    }

    vector<optional<double>> Founds::selong() const
    {
        int n = data.size();
        vector<optional<double>> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       {
                           Ephemeris eph = (found.ephemeris.num_vertices() == 1)
                                               ? found.ephemeris
                                               : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return eph.data(0).selong;
                       });
        return v;
    }

    vector<optional<double>> Founds::true_anomaly() const
    {
        int n = data.size();
        vector<optional<double>> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       {
                           Ephemeris eph = (found.ephemeris.num_vertices() == 1)
                                               ? found.ephemeris
                                               : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return eph.data(0).true_anomaly;
                       });
        return v;
    }

    vector<optional<double>> Founds::sangle() const
    {
        int n = data.size();
        vector<optional<double>> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       {
                           Ephemeris eph = (found.ephemeris.num_vertices() == 1)
                                               ? found.ephemeris
                                               : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return eph.data(0).sangle;
                       });
        return v;
    }

    vector<optional<double>> Founds::vangle() const
    {
        int n = data.size();
        vector<optional<double>> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       {
                           Ephemeris eph = (found.ephemeris.num_vertices() == 1)
                                               ? found.ephemeris
                                               : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return eph.data(0).vangle;
                       });
        return v;
    }

    vector<optional<double>> Founds::vmag() const
    {
        int n = data.size();
        vector<optional<double>> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       {
                           Ephemeris eph = (found.ephemeris.num_vertices() == 1)
                                               ? found.ephemeris
                                               : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return eph.data(0).vmag;
                       });
        return v;
    }

    json::array Founds::as_json()
    {
        json::array array;
        for (Found found : data)
            array.emplace_back(found.as_json());
        return array;
    }

    std::ostream &operator<<(std::ostream &os, const Found &found)
    {
        // if found.ephemeris is a segment, interpolate it to observation mid-time.
        Ephemeris eph;
        if (found.ephemeris.num_vertices() > 1)
        {
            double mjd = (found.observation.mjd_start() + found.observation.mjd_stop()) / 2;
            eph = found.ephemeris.interpolate(mjd);
        }
        else
            eph = found.ephemeris;

        os << found.observation << "  " << eph;

        return os;
    }

    std::ostream &operator<<(std::ostream &os, const Founds &founds)
    {
        Table table;
        table.add_column("observation_id", "%" PRId64, founds.observation_id());
        table.add_column("source", "%s", founds.source());
        table.add_column("product_id", "%s", founds.product_id());
        table.add_column("observatory", "%s", founds.observatory());
        table.add_column("mjd_start", "%.6lf", founds.mjd_start());
        table.add_column("mjd_stop", "%.6lf", founds.mjd_stop());
        table.add_column("exposure", "%.3lf", founds.exposure());
        if (founds.format.show_fov)
            table.add_column("fov", "%s", founds.fov());
        table.add_column("moving_target_id", "%" PRId64, founds.moving_target_id());
        table.add_column("designation", "%s", founds.designation());
        table.add_column("small_body", "%s", founds.small_body());
        table.add_column("mjd", "%.6lf", founds.mjd());
        table.add_column("tmtp", "%.6lf", founds.tmtp());
        table.add_column("ra", "%.6lf", founds.ra());
        table.add_column("dec", "%.6lf", founds.dec());
        table.add_column("mu", "%.2f", founds.mu());
        table.add_column("mu_theta", "%.3f", founds.mu_theta());
        table.add_column("rh", "%.4f", founds.rh());
        table.add_column("delta", "%.4f", founds.delta());
        table.add_column("phase", "%.3f", founds.phase());
        table.add_column("selong", "%.3f", founds.selong());
        table.add_column("true_anomaly", "%.3f", founds.true_anomaly());
        table.add_column("sangle", "%.3f", founds.sangle());
        table.add_column("vangle", "%.3f", founds.vangle());
        table.add_column("unc_a", "%.3f", founds.unc_a());
        table.add_column("unc_b", "%.3f", founds.unc_b());
        table.add_column("unc_th", "%.3f", founds.unc_theta());
        table.add_column("vmag", "%.3f", founds.vmag());

        os << table;
        return os;
    }
}