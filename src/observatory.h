#ifndef SBS_OBSERVATORY_H_
#define SBS_OBSERVATORY_H_

#include <cmath>
#include <string>
#include <s2/s2latlng.h>

#include "util.h"
#include "sofa/sofa.h"

using std::atan2;
using std::cos;
using std::sin;
using std::string;

namespace sbsearch
{
    struct Observatory
    {
        // Parallax constants: longitude, rho cos(phi), rho sin(phi)
        double longitude = 0; // deg E of Greenwich
        double rho_cos_phi = 0;
        double rho_sin_phi = 0;

        // Correct coordinates (RA, Dec), given date and distance (deg and au)
        // for parallax.
        //
        // Good to 10 arcsec at 1 Lunar distance.
        const S2LatLng parallax(const S2LatLng coords, const double mjd, const double delta) const
        {
            const double era = iauEra00(2400000.5, mjd);
            double ha = std::fmod(era + longitude * DEG - coords.lng().radians(), 2 * PI);
            if (ha < 0)
                ha += 2 * PI;

            // parallax
            // const double sin_pi = 4.26352098e-05 / delta; // sin(R_Earth / 1 au) / Delta
            const double sin_pi = 4.263521245426389e-05 / delta; // sin(R_Earth / 1 au) / Delta
            const double cos_dec = cos(coords.lat().radians());

            double y = -rho_cos_phi * sin_pi * sin(ha);
            const double x = cos_dec - rho_cos_phi * sin_pi * cos(ha);
            const double delta_ra = atan2(y, x);

            y = (sin(coords.lat().radians()) - rho_sin_phi * sin_pi) * cos(delta_ra);
            const double dec = atan2(y, x);

            return S2LatLng::FromRadians(dec, coords.lng().radians() + delta_ra);
        }
    };
}

#endif // SBS_OBSERVATORY_H_
