#ifndef SBS_EPHEMERIS_EXTRAPOLATE_H_
#define SBS_EPHEMERIS_EXTRAPOLATE_H_

#include "ephemeris/ephemeris.h"

namespace sbsearch::ephemeris
{
    // For ephemeris extrapolation: BACKWARDS to extrapolate before the
    // first vertex, FORWARDS to extrapolate beyond the last vertex.
    enum class Extrapolate : uint8
    {
        BACKWARDS,
        FORWARDS
    };

    // Linearly (on the sphere) extrapolate ephemeris by amount `distance`
    // in radians
    Ephemeris::Datum extrapolate(const Ephemeris::Data &data,
                                 const double distance,
                                 const Extrapolate direction);
    Ephemeris extrapolate(const Ephemeris &eph,
                          const double distance,
                          const Extrapolate direction);
}

#endif