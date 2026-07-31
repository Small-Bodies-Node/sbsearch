#ifndef SBS_EPHEMERIS_INTERPOLATE_H_
#define SBS_EPHEMERIS_INTERPOLATE_H_

#include "config.h"
#include "date.h"
#include "ephemeris.h"
#include "sbsearch.h"

using namespace sbsearch;

namespace sbsearch::sbs_ephemeris
{
    // helper for ephemeris interpolation: step is in units of days
    Ephemeris::Data interpolate(Ephemeris &eph, const Date start, const Date stop, const double step);
}

#endif