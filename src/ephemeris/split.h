#ifndef SBS_EPHEMERIS_SPLIT_H_
#define SBS_EPHEMERIS_SPLIT_H_

#include <vector>
#include "ephemeris/ephemeris.h"

using std::vector;

namespace sbsearch::ephemeris
{
    // Split ephemeris in segments of approximate length `length` in degrees and `time` in days.
    vector<Ephemeris> split(const Ephemeris &eph, const double length, const double time);
}

#endif
