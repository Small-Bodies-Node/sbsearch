#include "config.h"

#include "ephemeris.h"
#include "util/math.h"

namespace sbsearch
{
    Ephemeris::Datum Ephemeris::extrapolate(const double distance, Ephemeris::Extrapolate direction) const
    {
        int i, j;
        if (direction == Ephemeris::Extrapolate::BACKWARDS)
        {
            i = 1;
            j = 0;
        }
        else
        {
            i = num_vertices() - 2;
            j = num_vertices() - 1;
        }
        const Datum d1 = data_[i];
        const Datum d2 = data_[j];
        const S2Point p1 = d1.as_s2point();
        const S2Point p2 = d2.as_s2point();
        const double length = S1Angle(p1, p2).radians();
        const double frac = 1 + distance / length;

        S2Point extrapolated = S2::Interpolate(p1, p2, frac).Normalize();

        Datum d;
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
}