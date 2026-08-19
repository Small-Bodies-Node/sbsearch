#ifndef SBS_OBSERVATORY_H_
#define SBS_OBSERVATORY_H_

#include <map>
#include <string>
#include <tuple>
#include <s2/s2latlng.h>

using std::string;

namespace sbsearch
{
    // forward declaration for Observatories
    struct Observatory;

    typedef std::map<string, Observatory> Observatories;

    struct Observatory
    {
        // Parallax constants: longitude, rho cos(phi), rho sin(phi)
        double longitude = 0; // deg E of Greenwich
        double rho_cos_phi = 0;
        double rho_sin_phi = 0;
        string name = "";

        bool operator==(const Observatory &other) const
        {
            return (std::tie(longitude, rho_cos_phi, rho_sin_phi, name) ==
                    std::tie(other.longitude, other.rho_cos_phi, other.rho_sin_phi, other.name));
        }

        bool operator!=(const Observatory &other) const
        {
            return !(*this == other);
        }

        // Correct coordinates (RA, Dec), given date and distance (deg and au)
        // for parallax.
        //
        // Good to 10 arcsec at 1 Lunar distance.
        S2LatLng parallax(const S2LatLng &coords, const double mjd, const double delta) const;

        // Convenience method to insert this object into an Observatories
        // mapping.
        inline void insert_into(Observatories &observatories)
        {
            observatories[name] = *this;
        };
    };

}

#endif // SBS_OBSERVATORY_H_
