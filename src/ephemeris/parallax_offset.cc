#include "config.h"

#include "ephemeris/ephemeris.h"
#include "ephemeris/parallax_offset.h"
#include "observatory.h"

namespace sbsearch::ephemeris
{
    Ephemeris parallax_offset(const Ephemeris &eph, const Observatory &observatory)
    {
        Ephemeris::Data data(eph.data());
        for (int i = 0; i < data.size(); i++)
        {
            const S2LatLng coords = observatory.parallax(
                data[i].as_s2latlng(),
                data[i].mjd.value(),
                data[i].delta.value());
            data[i].radec(coords.Normalized().ToPoint());
        }
        return {eph.target(), data, eph.format};
    }
}