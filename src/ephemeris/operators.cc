#include "config.h"

#include <iostream>

#include "date.h"
#include "ephemeris.h"
#include "table.h"

using namespace sbsearch::table;

namespace sbsearch
{
    const Ephemeris Ephemeris::operator[](const int k) const
    {
        Ephemeris eph = Ephemeris(target_, {data(k)});
        return eph;
    }

    std::ostream &operator<<(std::ostream &os, const Ephemeris &ephemeris)
    {
        Table table;
        if (ephemeris.format.date == DateFormat::CALENDAR)
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

    bool
    Ephemeris::operator==(const Ephemeris &other) const
    {
        return ((target_ == other.target()) && (data_ == other.data()));
    }
}