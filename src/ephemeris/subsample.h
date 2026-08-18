#ifndef SBS_EPHEMERIS_SUBSAMPLE_H_
#define SBS_EPHEMERIS_SUBSAMPLE_H_

#include "ephemeris/ephemeris.h"

namespace sbsearch::ephemeris
{
    /* Get a section of the ephemeris based on the given date range.

    The ephemeris is interpolated to match the start and stop times.

    No data is returned for times not covered by the data.
    */
    Ephemeris::Data subsample(const Ephemeris::Data &data,
                              const double mjd_start,
                              const double mjd_stop);
    Ephemeris subsample(const Ephemeris &eph,
                        const double mjd_start,
                        const double mjd_stop);
}

#endif