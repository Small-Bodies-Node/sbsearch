#include "config.h"

#include "observatory.h"
#include "vertices.h"
#include "ephemeris/ephemeris.h"
#include "ephemeris/parallax_offset.h"

namespace sbsearch::ephemeris
{
    Ephemeris parallax_offset(const Ephemeris &eph, const Observatory &observatory)
    {
        Ephemeris::Data data(eph.data());
        for (auto &d : data)
        {
            const S2LatLng coords = observatory.parallax(
                make_latlng(d),
                d.mjd.value(),
                d.delta.value());
            d.radec(coords.Normalized().ToPoint());
        }
        return {eph.target(), data, eph.format};
    }
}