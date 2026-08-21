#include "config.h"

#include <cinttypes>
#include <optional>
#include <ostream>

#include "ephemeris/ephemeris.h"
#include "ephemeris/interpolate.h"
#include "found.h"
#include "observation.h"
#include "table.h"

using namespace sbsearch::table;

using sbsearch::ephemeris::Ephemeris;
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

    Founds::Founds(const vector<Found> &founds)
    {
        append(founds);
    }

    Founds::Founds(vector<Found> &&founds)
    {
        append(std::move(founds));
    }

    void Founds::append(const Found &found)
    {
        data.push_back(found);
    }

    void Founds::append(Found &&found)
    {
        data.push_back(std::move(found));
    }

    void Founds::append(const vector<Found> &founds)
    {
        data.reserve(data.size() + founds.size());
        std::copy(founds.begin(), founds.end(), std::back_inserter(data));
    }

    void Founds::append(vector<Found> &&founds)
    {
        data.reserve(data.size() + founds.size());
        std::move(founds.begin(), founds.end(), std::back_inserter(data));
    }

    void Founds::append(const Founds &founds)
    {
        append(founds.data);
    }

    void Founds::append(Founds &&founds)
    {
        append(std::move(founds.data));
    }

    vector<string> Founds::source() const
    {
        return get_observational_vector<vector<string>>(&Observation::source);
    };

    // Observatories.
    vector<string> Founds::observatory() const
    {
        return get_observational_vector<vector<string>>(&Observation::observatory);
    };

    // Observational product IDs.
    vector<string> Founds::product_id() const
    {
        return get_observational_vector<vector<string>>(&Observation::product_id);
    };

    // Observational fields-of-view.
    vector<string> Founds::fov() const
    {
        return get_observational_vector<vector<string>>(&Observation::fov);
    };

    // Observation IDs.
    vector<optional<int64_t>> Founds::observation_id() const
    {
        return get_observational_vector<vector<optional<int64_t>>>(&Observation::observation_id);
    };

    // Observation start dates.
    vector<double> Founds::mjd_start() const
    {
        return get_observational_vector<vector<double>>(&Observation::mjd_start);
    };

    // Observation stop dates.
    vector<double> Founds::mjd_stop() const
    {
        return get_observational_vector<vector<double>>(&Observation::mjd_stop);
    };

    // Observation exposure times (s).
    vector<double> Founds::exposure() const
    {
        return get_observational_vector<vector<double>>(&Observation::exposure);
    };

    // Found moving target IDs.
    vector<optional<int64_t>> Founds::moving_target_id() const
    {
        return get_target_vector<vector<optional<int64_t>>>(&MovingTarget::moving_target_id);
    };

    // Found moving target designations.
    vector<string> Founds::designation() const
    {
        return get_target_vector<vector<string>>(&MovingTarget::designation);
    };

    // Found moving target small-body flags.
    vector<bool> Founds::small_body() const
    {
        return get_target_vector<vector<bool>>(&MovingTarget::small_body);
    };

    // Times at which the ephemerides are calculated.
    vector<double> Founds::mjd() const
    {
        return get_ephemeral_vector<vector<double>>(&Ephemeris::Datum::mjd);
    };

    vector<string> Founds::date() const
    {
        vector<double> mjds = mjd();
        vector<string> v(size());
        std::transform(mjds.begin(), mjds.end(), v.begin(), Date::MJDToCalendar);
        return v;
    }
    // Approximate time relative to perihelion (days).
    vector<optional<double>> Founds::tmtp() const
    {
        return get_ephemeral_vector<vector<optional<double>>>(&Ephemeris::Datum::tmtp);
    };

    // Ephemeris right ascension.
    vector<double> Founds::ra() const
    {
        return get_ephemeral_vector<vector<double>>(&Ephemeris::Datum::ra);
    };

    // Ephemeris declination.
    vector<double> Founds::dec() const
    {
        return get_ephemeral_vector<vector<double>>(&Ephemeris::Datum::dec);
    };

    // Ephemeris uncertainty ellipse semi-major axis.
    vector<optional<double>> Founds::unc_a() const
    {
        return get_ephemeral_vector<vector<optional<double>>>(&Ephemeris::Datum::unc_a);
    };

    // Ephemeris uncertainty ellipse semi-minor axis.
    vector<optional<double>> Founds::unc_b() const
    {
        return get_ephemeral_vector<vector<optional<double>>>(&Ephemeris::Datum::unc_b);
    };

    // Ephemeris proper motion.
    vector<double> Founds::mu() const
    {
        return get_ephemeral_vector<vector<double>>(&Ephemeris::Datum::mu);
    };

    // Ephemeris proper motion direction.
    vector<double> Founds::mu_theta() const
    {
        return get_ephemeral_vector<vector<double>>(&Ephemeris::Datum::mu_theta);
    };

    // Ephemeris uncertainty ellipse semi-major axis position angle (deg E of N).
    vector<optional<double>> Founds::unc_theta() const
    {
        return get_ephemeral_vector<vector<optional<double>>>(&Ephemeris::Datum::unc_theta);
    };

    // Heliocentric distance (au).
    vector<double> Founds::rh() const
    {
        return get_ephemeral_vector<vector<double>>(&Ephemeris::Datum::rh);
    };

    // Observer-target distance (au).
    vector<double> Founds::delta() const
    {
        return get_ephemeral_vector<vector<double>>(&Ephemeris::Datum::delta);
    };

    // Sun-target-observer angle (deg).
    vector<double> Founds::phase() const
    {
        return get_ephemeral_vector<vector<double>>(&Ephemeris::Datum::phase);
    };

    // Solar elongation (deg).
    vector<optional<double>> Founds::selong() const
    {
        return get_ephemeral_vector<vector<optional<double>>>(&Ephemeris::Datum::selong);
    };

    // Orbital true anomaly (deg).
    vector<optional<double>> Founds::true_anomaly() const
    {
        return get_ephemeral_vector<vector<optional<double>>>(&Ephemeris::Datum::true_anomaly);
    };

    // Projected target-sun vector position angle (deg).
    vector<optional<double>> Founds::sangle() const
    {
        return get_ephemeral_vector<vector<optional<double>>>(&Ephemeris::Datum::sangle);
    };

    // Projected target velocity vector position angle (deg).
    vector<optional<double>> Founds::vangle() const
    {
        return get_ephemeral_vector<vector<optional<double>>>(&Ephemeris::Datum::vangle);
    };

    // Target apparent magnitude (meaning varies depending on ephemeris source).
    vector<optional<double>> Founds::vmag() const
    {
        return get_ephemeral_vector<vector<optional<double>>>(&Ephemeris::Datum::vmag);
    };

    vector<double> Founds::mjd_added() const
    {
        int n = data.size();
        vector<double> v(n);

        auto get = std::mem_fn(&Found::mjd_added);
        std::transform(data.begin(), data.end(), v.begin(), get);
        return v;
    }

    std::ostream &operator<<(std::ostream &os, const Found &found)
    {
        // if found.ephemeris is a segment, interpolate it to observation mid-time.
        Ephemeris::Datum point;
        if (found.ephemeris.num_vertices() > 1)
        {
            double mjd = (found.observation.mjd_start() + found.observation.mjd_stop()) / 2;
            point = ephemeris::interpolate(found.ephemeris.data(), mjd);
        }
        else
            point = found.ephemeris.data(0);

        os << found.observation << "  " << Ephemeris(found.ephemeris.target(), {point});

        return os;
    }

    std::ostream &operator<<(std::ostream &os, const Founds &founds)
    {
        auto format_dates = [](vector<double> mjds)
        {
            vector<string> dates(mjds.size());
            std::transform(mjds.begin(), mjds.end(), dates.begin(),
                           [](const double mjd)
                           { return Date(mjd).iso(); });
            return dates;
        };

        Table table;
        table.add(Column("observation_id", "%" PRId64, founds.observation_id()));
        table.add(Column("source", "%s", founds.source()));
        table.add(Column("product_id", "%s", founds.product_id()));
        table.add(Column("observatory", "%s", founds.observatory()));

        if (founds.observation_format.date == Date::Format::MJD)
        {
            table.add(Column("date_start", "%.6lf", founds.mjd_start()));
            table.add(Column("date_stop", "%.6lf", founds.mjd_stop()));
        }
        else
        {
            table.add(Column("date_start", "%s", format_dates(founds.mjd_start())));
            table.add(Column("date_stop", "%s", format_dates(founds.mjd_stop())));
        }

        table.add(Column("exposure", "%.3lf", founds.exposure()));

        if (founds.observation_format.show_fov)
            table.add(Column("fov", "%s", founds.fov()));

        table.add(Column("moving_target_id", "%" PRId64, founds.moving_target_id()));
        table.add(Column("designation", "%s", founds.designation()));
        table.add(Column("small_body", "%s", founds.small_body()));

        if (founds.ephemeris_format.date == Date::Format::MJD)
            table.add(Column("date", "%.6lf", founds.mjd()));
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