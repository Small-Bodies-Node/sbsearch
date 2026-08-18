#include <algorithm>

#include "config.h"
#include "ephemeris/ephemeris.h"
#include "ephemeris/interpolate.h"
#include "ephemeris/subsample.h"

namespace sbsearch::ephemeris
{
    Ephemeris::Data subsample(const Ephemeris::Data &data, const double mjd_start, const double mjd_stop)
    {
        // Is there data to subsample>?  Is the requested time period outside of the data?
        if (data.size() == 0 || (data.front().mjd > mjd_stop) || (data.back().mjd < mjd_start))
            return {};

        Ephemeris::Data result;

        // find any whole segments between start and end
        auto end = data.cend();
        auto next = std::lower_bound(data.cbegin(), end, mjd_start,
                                     [](const Ephemeris::Datum &d, const double &mjd_start)
                                     { return d.mjd.value() < mjd_start; });
        auto last = std::upper_bound(next, end, mjd_stop,
                                     [](const double &mjd_stop, const Ephemeris::Datum &d)
                                     { return mjd_stop <= d.mjd.value(); });

        result.reserve(last - next + 2);

        // Was interpolation between two epochs requested?
        if (next->mjd > mjd_start)
            result.push_back(interpolate(data, mjd_start));

        // Are there points between next and last?
        if (next != end && next != last)
            std::copy(next, last, std::back_inserter(result));

        // Is last the end or was interpolation between two epochs requested?
        if (last->mjd == mjd_stop)
            result.push_back(*last);
        else
            result.push_back(interpolate(data, mjd_stop));

        return result;
    }

    Ephemeris subsample(const Ephemeris &eph, const double mjd_start, const double mjd_stop)
    {
        return {eph.target(), subsample(eph.data(), mjd_start, mjd_stop), eph.format};
    }
}