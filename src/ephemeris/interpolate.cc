#include <algorithm>
#include <functional>
#include <valarray>
#include <vector>

#include "config.h"
#include "vertices.h"
#include "ephemeris/ephemeris.h"
#include "ephemeris/interpolate.h"
#include "util/math.h"
#include "util/optional.h"

using std::cerr;
using std::endl;
using std::valarray;
using std::vector;

namespace sbsearch::ephemeris
{
    Ephemeris::Datum interpolate(const Ephemeris::Data &data, const double mjd0)
    {
        if ((mjd0 < data.front().mjd) || (mjd0 > data.back().mjd))
            throw std::runtime_error("Interpolation beyond ephemeris time range not supported");

        // order of polynomial interpolation is n - 1
        const gsl_interp_type *interp_type;
        int n = std::min(4, (int)data.size());
        if (n < 2)
            throw std::runtime_error("Cannot interpolate with an ephemeris of length " + std::to_string(n));
        else if (n == 2)
            interp_type = gsl_interp_linear;
        else
            interp_type = gsl_interp_polynomial;

        // find the nearest segment; left and right are inclusive
        auto right = std::lower_bound(data.cbegin(), data.cend(), mjd0,
                                      [](const Ephemeris::Datum &d, const double mjd)
                                      { return d.mjd < mjd; });
        auto left = right - n + 1;

        // shift to center on mjd, taking care with integer math
        std::advance(right, (n - 2) / 2);
        std::advance(left, (n - 2) / 2);

        if (left < data.cbegin())
        {
            int offset = std::distance(left, data.cbegin());
            std::advance(right, offset);
            std::advance(left, offset);
        }

        if (right >= data.cend())
        {
            int offset = std::distance(right, data.cend() - 1);
            std::advance(right, offset);
            std::advance(left, offset);
        }

        if ((left < data.cbegin()) || right >= data.cend())
            throw std::runtime_error("Interpolation range is out of bounds.");

        // use the limits to define the segment that we will interpolate
        const Ephemeris segment({}, {left, right + 1});
        const vector<double> mjd = segment.mjd();

        unique_interp_accel_ptr accel(gsl_interp_accel_alloc(), &gsl_interp_accel_free);
        unique_interp_ptr interp(gsl_interp_alloc(interp_type, n), &gsl_interp_free);

        Ephemeris::Datum d;
        d.mjd = mjd0;
        d.tmtp = interpolate_optional_vector_(mjd0, mjd, segment.tmtp(), interp.get(), accel.get());
        d.mu = interpolate_vector_(mjd0, mjd, segment.mu(), interp.get(), accel.get());
        d.mu_theta = interpolate_vector_(mjd0, mjd, segment.mu_theta(), interp.get(), accel.get(), true);
        d.unc_a = interpolate_optional_vector_(mjd0, mjd, segment.unc_a(), interp.get(), accel.get());
        d.unc_b = interpolate_optional_vector_(mjd0, mjd, segment.unc_b(), interp.get(), accel.get());
        d.unc_theta = interpolate_optional_vector_(mjd0, mjd, segment.unc_theta(), interp.get(), accel.get(), true);
        d.rh = interpolate_vector_(mjd0, mjd, segment.rh(), interp.get(), accel.get());
        d.delta = interpolate_vector_(mjd0, mjd, segment.delta(), interp.get(), accel.get());
        d.phase = interpolate_vector_(mjd0, mjd, segment.phase(), interp.get(), accel.get(), true);
        d.selong = interpolate_optional_vector_(mjd0, mjd, segment.selong(), interp.get(), accel.get(), true);
        d.true_anomaly = interpolate_optional_vector_(mjd0, mjd, segment.true_anomaly(), interp.get(), accel.get(), true);
        d.sangle = interpolate_optional_vector_(mjd0, mjd, segment.sangle(), interp.get(), accel.get(), true);
        d.vangle = interpolate_optional_vector_(mjd0, mjd, segment.vangle(), interp.get(), accel.get(), true);
        d.vmag = interpolate_optional_vector_(mjd0, mjd, segment.vmag(), interp.get(), accel.get());

        // for ra and dec, convert to Cartesian coordinates, then interpolate, and convert back
        vector<double> vx, vy, vz;
        for (auto const &p : segment.data())
        {
            // in limited testing, multiplying by delta did not significantly
            // improve things
            auto point = make_point(p);
            vx.push_back(point.x());
            vy.push_back(point.y());
            vz.push_back(point.z());
        }

        double x = interpolate_vector_(mjd0, mjd, vx, interp.get(), accel.get());
        double y = interpolate_vector_(mjd0, mjd, vy, interp.get(), accel.get());
        double z = interpolate_vector_(mjd0, mjd, vz, interp.get(), accel.get());

        S2LatLng point(S2Point(x, y, z).Normalize());

        d.ra = point.lng().degrees();
        d.dec = point.lat().degrees();

        return d;
    }

    Ephemeris interpolate(const Ephemeris &eph, const double mjd0)
    {
        return {eph.target(), {interpolate(eph.data(), mjd0)}, eph.format};
    }

    double interpolate_valarray_(const double x0,
                                 const vector<double> &x,
                                 const valarray<double> &y,
                                 gsl_interp *interp,
                                 gsl_interp_accel *accel,
                                 bool angle)
    {
        double offset = 0;
        valarray<double> yy(y);
        if (angle)
        {
            offset = y.sum() / y.size();
            yy -= offset;
        }

        gsl_interp_init(interp, x.data(), &yy[0], yy.size());
        double y0 = gsl_interp_eval(interp, x.data(), &yy[0], x0, accel);
        return y0 + offset;
    }

    double interpolate_vector_(const double x0,
                               const vector<double> &x,
                               const vector<double> &y,
                               gsl_interp *interp,
                               gsl_interp_accel *accel,
                               bool angle)
    {
        return interpolate_valarray_(x0, x, {y.data(), y.size()}, interp, accel, angle);
    }

    optional<double> interpolate_optional_vector_(const double x0,
                                                  const vector<double> &x,
                                                  const vector<optional<double>> &yo,
                                                  gsl_interp *interp,
                                                  gsl_interp_accel *accel,
                                                  bool angle)
    {
        optional<double> result;

        try
        {
            vector<double> y = util::optionals_to_values(yo);
            result = interpolate_valarray_(x0, x, {y.data(), y.size()}, interp, accel, angle);
        }
        catch (std::bad_optional_access)
        {
            result = std::nullopt;
        }

        return result;
    }
}