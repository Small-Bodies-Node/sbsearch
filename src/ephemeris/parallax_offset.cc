#include "config.h"

#include "ephemeris.h"

namespace sbsearch
{
    Ephemeris Ephemeris::parallax_offset(const Observatory &observatory)
    {
        Data new_data(data_);
        for (int i = 0; i < data_.size(); i++)
        {
            const S2LatLng coords = observatory.parallax(
                new_data[i].as_s2latlng(),
                new_data[i].mjd.value(),
                new_data[i].delta.value());
            new_data[i].radec(coords.Normalized().ToPoint());
        }
        return Ephemeris(target_, new_data);
    }
}