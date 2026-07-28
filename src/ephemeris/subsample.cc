#include "config.h"

#include "ephemeris.h"

namespace sbsearch
{
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
}