#ifndef SBS_POLYLINE_H_
#define SBS_POLYLINE_H_

#include <s2/s2polyline.h>
#include "ephemeris/ephemeris.h"

using sbsearch::ephemeris::Ephemeris;

namespace sbsearch
{
    // Represent the ephemeris as a polyline.
    inline S2Polyline make_polyline(const Ephemeris eph)
    {
        return S2Polyline(eph.vertices());
    }
}

#endif
