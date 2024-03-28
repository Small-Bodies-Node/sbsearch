#include "config.h"

#include <ostream>

#include "ephemeris.h"
#include "found.h"
#include "observation.h"
#include "table.h"

using sbsearch::table::Table;

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

    vector<int64> Founds::observation_id() const
    {
        int n = data.size();
        vector<int64> v(n);
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

    vector<int64> Founds::moving_target_id() const
    {
        int n = data.size();
        vector<int64> v(n);
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

    vector<double> Founds::mjd() const
    {
        int n = data.size();
        vector<double> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       {
                           Ephemeris eph = (found.ephemeris.num_vertices() == 1)
                                               ? found.ephemeris
                                               : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return eph.data(0).mjd;
                       });
        return v;
    }

    vector<double> Founds::tmtp() const
    {
        int n = data.size();
        vector<double> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       {
                           Ephemeris eph = (found.ephemeris.num_vertices() == 1)
                                               ? found.ephemeris
                                               : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return eph.data(0).tmtp;
                       });
        return v;
    }

    vector<double> Founds::ra() const
    {
        int n = data.size();
        vector<double> v(n);
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

    vector<double> Founds::dec() const
    {
        int n = data.size();
        vector<double> v(n);
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

    vector<double> Founds::unc_a() const
    {
        int n = data.size();
        vector<double> v(n);
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

    vector<double> Founds::unc_b() const
    {
        int n = data.size();
        vector<double> v(n);
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

    vector<double> Founds::unc_theta() const
    {
        int n = data.size();
        vector<double> v(n);
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

    vector<double> Founds::rh() const
    {
        int n = data.size();
        vector<double> v(n);
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

    vector<double> Founds::delta() const
    {
        int n = data.size();
        vector<double> v(n);
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

    vector<double> Founds::phase() const
    {
        int n = data.size();
        vector<double> v(n);
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

    vector<double> Founds::selong() const
    {
        int n = data.size();
        vector<double> v(n);
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

    vector<double> Founds::true_anomaly() const
    {
        int n = data.size();
        vector<double> v(n);
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

    vector<double> Founds::sangle() const
    {
        int n = data.size();
        vector<double> v(n);
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

    vector<double> Founds::vangle() const
    {
        int n = data.size();
        vector<double> v(n);
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

    vector<double> Founds::vmag() const
    {
        int n = data.size();
        vector<double> v(n);
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
        bool show_fov = std::max_element(founds.begin(), founds.end(),
                                         [](const Found &a, const Found &b)
                                         { return a.observation.format.show_fov < b.observation.format.show_fov; })
                            ->observation.format.show_fov;

        Table table;
        table.add_column("observation_id", "%" PRId64, founds.observation_id());
        table.add_column("source", "%s", founds.source());
        table.add_column("observatory", "%s", founds.observatory());
        table.add_column("mjd_start", "%.6lf", founds.mjd_start());
        table.add_column("mjd_stop", "%.6lf", founds.mjd_stop());
        table.add_column("expsoure", "%.3lf", founds.exposure());
        if (show_fov)
            table.add_column("fov", "%s", founds.fov());
        table.add_column("moving_target_id", "%" PRId64, founds.moving_target_id());
        table.add_column("designation", "%s", founds.designation());
        table.add_column("small_body", "%s", founds.small_body());
        table.add_column("mjd", "%.6lf", founds.mjd());
        table.add_column("tmtp", "%.6lf", founds.tmtp());
        table.add_column("ra", "%.6lf", founds.ra());
        table.add_column("dec", "%.6lf", founds.dec());
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