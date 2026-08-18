#ifndef SBS_SBSEPHEMERIS_INTERPOLATE_H_
#define SBS_SBSEPHEMERIS_INTERPOLATE_H_

#include "config.h"
#include "date.h"
#include "ephemeris/ephemeris.h"
#include "sbsearch.h"

using namespace sbsearch;
using sbsearch::ephemeris::Ephemeris;

namespace sbsearch::sbs_ephemeris
{
    // helper for ephemeris interpolation: step is in units of days
    Ephemeris::Data interpolate(Ephemeris &eph, const Date start, const Date stop, const double step);
}

#endif