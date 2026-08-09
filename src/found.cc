#include "config.h"

#include <cinttypes>
#include <optional>
#include <ostream>
#include <boost/json.hpp>

#include "ephemeris.h"
#include "found.h"
#include "observation.h"
#include "table.h"

using namespace sbsearch::table;
using std::optional;

namespace sbsearch
{
    const bool Found::operator==(const Found &other) const
    {
        return ((observation == other.observation) && (ephemeris == other.ephemeris));
    };

    const bool Found::operator!=(const Found &other) const
    {
        return !(*this == other);
    };

    json::object Found::as_json()
    {
        json::object obj;

        // if found.ephemeris is a segment, interpolate it to observation mid-time.
        Ephemeris eph(ephemeris.target(), {});
        if (ephemeris.num_vertices() > 1)
        {
            double mjd = (observation.mjd_start() + observation.mjd_stop()) / 2;
            eph.append(ephemeris.interpolate(mjd));
        }
        else
            eph = ephemeris;

        for (auto item : observation.as_json())
            obj[item.key()] = item.value();

        json::object eph_object = *eph.as_json().at(0).if_object();
        for (auto item : eph_object)
            obj[item.key()] = item.value();

        obj["mjd_added"] = mjd_added;

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
                           auto point = (found.ephemeris.num_vertices() == 1)
                                            ? found.ephemeris.data(0)
                                            : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return point.mjd.value_or(-1);
                       });
        return v;
    }

    vector<string> Founds::date() const
    {
        int n = data.size();
        vector<string> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       {
                           auto point = (found.ephemeris.num_vertices() == 1)
                                            ? found.ephemeris.data(0)
                                            : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return point.mjd ? Date(point.mjd.value()).iso() : "null";
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
                           auto point = (found.ephemeris.num_vertices() == 1)
                                            ? found.ephemeris.data(0)
                                            : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return point.tmtp.value_or(-1);
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
                           auto point = (found.ephemeris.num_vertices() == 1)
                                            ? found.ephemeris.data(0)
                                            : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return point.ra;
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
                           auto point = (found.ephemeris.num_vertices() == 1)
                                            ? found.ephemeris.data(0)
                                            : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return point.dec;
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
                           auto point = (found.ephemeris.num_vertices() == 1)
                                            ? found.ephemeris.data(0)
                                            : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return point.mu;
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
                           auto point = (found.ephemeris.num_vertices() == 1)
                                            ? found.ephemeris.data(0)
                                            : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return point.mu_theta;
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
                           auto point = (found.ephemeris.num_vertices() == 1)
                                            ? found.ephemeris.data(0)
                                            : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return point.unc_a;
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
                           auto point = (found.ephemeris.num_vertices() == 1)
                                            ? found.ephemeris.data(0)
                                            : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return point.unc_b;
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
                           auto point = (found.ephemeris.num_vertices() == 1)
                                            ? found.ephemeris.data(0)
                                            : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return point.unc_theta;
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
                           auto point = (found.ephemeris.num_vertices() == 1)
                                            ? found.ephemeris.data(0)
                                            : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return point.rh;
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
                           auto point = (found.ephemeris.num_vertices() == 1)
                                            ? found.ephemeris.data(0)
                                            : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return point.delta;
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
                           auto point = (found.ephemeris.num_vertices() == 1)
                                            ? found.ephemeris.data(0)
                                            : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return point.phase;
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
                           auto point = (found.ephemeris.num_vertices() == 1)
                                            ? found.ephemeris.data(0)
                                            : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return point.selong;
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
                           auto point = (found.ephemeris.num_vertices() == 1)
                                            ? found.ephemeris.data(0)
                                            : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return point.true_anomaly;
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
                           auto point = (found.ephemeris.num_vertices() == 1)
                                            ? found.ephemeris.data(0)
                                            : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return point.sangle;
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
                           auto point = (found.ephemeris.num_vertices() == 1)
                                            ? found.ephemeris.data(0)
                                            : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return point.vangle;
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
                           auto point = (found.ephemeris.num_vertices() == 1)
                                            ? found.ephemeris.data(0)
                                            : found.ephemeris.interpolate(found.observation.mjd_mid());
                           return point.vmag;
                       });
        return v;
    }

    vector<double> Founds::mjd_added() const
    {
        int n = data.size();
        vector<double> v(n);
        std::transform(data.begin(), data.end(), v.begin(),
                       [](const Found &found)
                       { return found.mjd_added; });
        return v;
    }
    json::array Founds::as_json() const
    {
        json::array array;
        for (Found found : data)
            array.emplace_back(found.as_json());
        return array;
    }

    std::ostream &operator<<(std::ostream &os, const Found &found)
    {
        // if found.ephemeris is a segment, interpolate it to observation mid-time.
        Ephemeris::Datum point;
        if (found.ephemeris.num_vertices() > 1)
        {
            double mjd = (found.observation.mjd_start() + found.observation.mjd_stop()) / 2;
            point = found.ephemeris.interpolate(mjd);
        }
        else
            point = found.ephemeris.data(0);

        os << found.observation << "  " << Ephemeris(found.ephemeris.target(), {point});

        return os;
    }

    std::ostream &operator<<(std::ostream &os, const Founds &founds)
    {
        Table table;
        table.add(Column("observation_id", "%" PRId64, founds.observation_id()));
        table.add(Column("source", "%s", founds.source()));
        table.add(Column("product_id", "%s", founds.product_id()));
        table.add(Column("observatory", "%s", founds.observatory()));
        table.add(Column("mjd_start", "%.6lf", founds.mjd_start()));
        table.add(Column("mjd_stop", "%.6lf", founds.mjd_stop()));
        table.add(Column("exposure", "%.3lf", founds.exposure()));
        if (founds.observation_format.show_fov)
            table.add(Column("fov", "%s", founds.fov()));
        table.add(Column("moving_target_id", "%" PRId64, founds.moving_target_id()));
        table.add(Column("designation", "%s", founds.designation()));
        table.add(Column("small_body", "%s", founds.small_body()));
        if (founds.ephemeris_format.date == Date::Format::MJD)
            table.add(Column("mjd", "%.6lf", founds.mjd()));
        else
            table.add(Column("date", "%19s", founds.date()));
        table.add(Column("tmtp", "%.6lf", founds.tmtp()));
        table.add(Column("ra", "%.6lf", founds.ra()));
        table.add(Column("dec", "%.6lf", founds.dec()));
        table.add(Column("mu", "%.2f", founds.mu()));
        table.add(Column("mu_theta", "%.3f", founds.mu_theta()));
        table.add(Column("rh", "%.4f", founds.rh()));
        table.add(Column("delta", "%.4f", founds.delta()));
        table.add(Column("phase", "%.3f", founds.phase()));
        table.add(Column("selong", "%.3f", founds.selong()));
        table.add(Column("true_anomaly", "%.3f", founds.true_anomaly()));
        table.add(Column("sangle", "%.3f", founds.sangle()));
        table.add(Column("vangle", "%.3f", founds.vangle()));
        table.add(Column("unc_a", "%.3f", founds.unc_a()));
        table.add(Column("unc_b", "%.3f", founds.unc_b()));
        table.add(Column("unc_th", "%.3f", founds.unc_theta()));
        table.add(Column("vmag", "%.3f", founds.vmag()));

        os << table;
        return os;
    }
}