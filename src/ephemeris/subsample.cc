#include "config.h"

#include "ephemeris.h"

namespace sbsearch
{
    Ephemeris Ephemeris::subsample(const double mjd_start, const double mjd_stop) const
    {
        Ephemeris eph(target_, {});

        // find any whole segments between start and end
        auto end = data_.cend();
        auto next = std::lower_bound(data_.cbegin(), end, mjd_start,
                                     [](const Datum &d, const double &mjd)
                                     { return d.mjd.value() < mjd; });
        auto last = std::upper_bound(next, end, mjd_stop,
                                     [](const double &mjd, const Datum &d)
                                     { return mjd < d.mjd.value(); });

        // Was interpolation between two epochs requested?
        if ((*next).mjd > mjd_start)
            eph.append(interpolate(mjd_start));

        // Are there points between the start and end?
        if (next != end && next != last)
            eph.append({next, last});

        // Was interpolation between two epochs requested?
        if (eph.data().back().mjd < mjd_stop)
            eph.append(interpolate(mjd_stop));

        return eph;
    }
}