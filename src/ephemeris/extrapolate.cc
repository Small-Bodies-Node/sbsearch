#include <s2/s2edge_distances.h>

#include "config.h"
#include "vertices.h"
#include "ephemeris/ephemeris.h"
#include "ephemeris/extrapolate.h"
#include "util/math.h"

namespace sbsearch::ephemeris
{
    Ephemeris::Datum extrapolate(const Ephemeris::Data &data,
                                 const double distance,
                                 const Extrapolate direction)
    {
        int i, j;
        if (direction == Extrapolate::BACKWARDS)
        {
            i = 1;
            j = 0;
        }
        else
        {
            i = data.size() - 2;
            j = data.size() - 1;
        }
        const Ephemeris::Datum d1 = data[i];
        const Ephemeris::Datum d2 = data[j];
        const S2Point p1 = make_point(d1);
        const S2Point p2 = make_point(d2);
        const double length = S1Angle(p1, p2).radians();
        const double frac = 1 + distance / length;

        S2Point extrapolated = S2::Interpolate(p1, p2, frac).Normalize();

        Ephemeris::Datum d;
        d.mjd = util::interp(d1.mjd, d2.mjd, frac);
        d.tmtp = util::interp(d1.tmtp, d2.tmtp, frac);
        d.radec(extrapolated);
        d.mu = util::interp(d1.mu, d2.mu, frac);
        d.mu_theta = util::interp(d1.mu_theta, d2.mu_theta, frac);
        d.unc_a = util::interp(d1.unc_a, d2.unc_a, frac);
        d.unc_b = util::interp(d1.unc_b, d2.unc_b, frac);
        d.unc_theta = util::interp(d1.unc_theta, d2.unc_theta, frac);
        d.rh = util::interp(d1.rh, d2.rh, frac);
        d.delta = util::interp(d1.delta, d2.delta, frac);
        d.phase = util::interp(d1.phase, d2.phase, frac);
        d.selong = util::interp(d1.selong, d2.selong, frac);
        d.true_anomaly = util::interp(d1.true_anomaly, d2.true_anomaly, frac);
        d.sangle = util::interp(d1.sangle, d2.sangle, frac);
        d.vangle = util::interp(d1.vangle, d2.vangle, frac);
        d.vmag = util::interp(d1.vmag, d2.vmag, frac);

        return d;
    }

    Ephemeris extrapolate(const Ephemeris &eph,
                          const double distance,
                          const Extrapolate direction)
    {
        return {eph.target(), {extrapolate(eph.data(), distance, direction)}, eph.format};
    }
}