#ifndef SBS_EPHEMERIS_SAMPLER_H_
#define SBS_EPHEMERIS_SAMPLER_H_

#include <shared_mutex>
#include "moving_target.h"
#include "ephemeris/ephemeris.h"

namespace sbsearch::ephemeris
{
    // Thread-safe ephemeris subsampling used for observation testing.
    class SafeSampler
    {
    public:
        // Append data onto the ephemeris.
        //
        // If the target is not defined, then it will be automatically updated.
        //
        // Unlike Ephemeris::append, if the last point of the current ephemeris
        // and the first point of the new ephemeris are the same, then the first
        // point of the ephemeris is skipped.
        void append(const Ephemeris &eph);

        // Subsample the ephemeris.
        Ephemeris subsample(double mjd_start, double mjd_stop) const;

    private:
        Ephemeris eph_;
        mutable std::shared_mutex access;
    };
}

#endif