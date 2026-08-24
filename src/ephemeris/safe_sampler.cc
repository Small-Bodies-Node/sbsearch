#include <mutex>
#include <shared_mutex>

#include "ephemeris/safe_sampler.h"
#include "ephemeris/ephemeris.h"
#include "ephemeris/subsample.h"

namespace sbsearch::ephemeris
{
    void SafeSampler::append(const Ephemeris &eph)
    {
        std::unique_lock<std::shared_mutex> lock(access);

        // Update moving target if the ID is not defined but the new ID is defined.
        if (!eph_.target().moving_target_id().has_value() && eph.target().moving_target_id().has_value())
            eph_.target(eph.target());

        if ((eph_.num_vertices() > 0) && (eph_.data(-1).mjd <= eph.data(0).mjd))
            eph_.append(eph.slice(1));
        else
            eph_.append(eph);
    }

    Ephemeris SafeSampler::subsample(double mjd_start, double mjd_stop) const
    {
        std::shared_lock<std::shared_mutex> lock(access);
        return ephemeris::subsample(eph_, mjd_start, mjd_stop);
    }
}