#ifndef SBS_PARALLAX_OFFSET_H_
#define SBS_PARALLAX_OFFSET_H_

#include "ephemeris/ephemeris.h"
#include "observatory.h"

namespace sbsearch::ephemeris
{
    // Offset the ephemeris for parallax.
    Ephemeris parallax_offset(const Ephemeris &eph, const Observatory &observatory);
}

#endif