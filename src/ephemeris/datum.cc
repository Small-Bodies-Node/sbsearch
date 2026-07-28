#include "config.h"

#include <tuple>
#include <boost/json.hpp>

#include "ephemeris.h"

namespace json = boost::json;

namespace sbsearch
{
    bool Ephemeris::Datum::operator==(const Ephemeris::Datum &other) const
    {
        return (std::tie(mjd, tmtp, ra, dec, mu, mu_theta,
                         unc_a, unc_b, unc_theta,
                         rh, delta, phase, selong, true_anomaly,
                         sangle, vangle, vmag) ==
                std::tie(other.mjd, other.tmtp, other.ra, other.dec, other.mu, other.mu_theta,
                         other.unc_a, other.unc_b, other.unc_theta,
                         other.rh, other.delta, other.phase, other.selong, other.true_anomaly,
                         other.sangle, other.vangle, other.vmag));
    }

    bool Ephemeris::Datum::operator!=(const Ephemeris::Datum &other) const
    {
        return !(*this == other);
    }

    json::object Ephemeris::Datum::as_json()
    {
        json::object datum;
        datum["mjd"] = json::value_from(mjd);
        datum["tmtp"] = json::value_from(tmtp);
        datum["ra"] = json::value_from(ra);
        datum["dec"] = json::value_from(dec);
        datum["mu"] = json::value_from(mu);
        datum["mu_theta"] = json::value_from(mu_theta);
        datum["unc_a"] = json::value_from(unc_a);
        datum["unc_b"] = json::value_from(unc_b);
        datum["unc_theta"] = json::value_from(unc_theta);
        datum["rh"] = json::value_from(rh);
        datum["delta"] = json::value_from(delta);
        datum["phase"] = json::value_from(phase);
        datum["selong"] = json::value_from(selong);
        datum["true_anomaly"] = json::value_from(true_anomaly);
        datum["sangle"] = json::value_from(sangle);
        datum["vangle"] = json::value_from(vangle);
        datum["vmag"] = json::value_from(vmag);
        return datum;
    }
}