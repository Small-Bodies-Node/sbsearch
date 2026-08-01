#include <algorithm>
#include <functional>
#include <valarray>
#include <vector>

#include "config.h"
#include "ephemeris.h"
#include "util/math.h"
#include "util/optional.h"

using std::cerr;
using std::endl;
using std::valarray;
using std::vector;

// Support RAII with GNU Scientific Library pointers
using unique_interp_accel_ptr = std::unique_ptr<gsl_interp_accel, decltype(&gsl_interp_accel_free)>;
using unique_interp_ptr = std::unique_ptr<gsl_interp, decltype(&gsl_interp_free)>;

namespace sbsearch
{
    // Interpolate a vector as a function of x to a specific date.
    double interpolate_vector_(const double x0,
                               const vector<double> &x,
                               const valarray<double> &y,
                               gsl_interp *interp,
                               gsl_interp_accel *accel,
                               bool angle = false)
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

    // Interpoalte a vector of optional values, and return nullopt if any values
    // are not defined.
    optional<double> interpolate_optional_vector_(const double x0,
                                                  const vector<double> &x,
                                                  const vector<optional<double>> &yo,
                                                  gsl_interp *interp,
                                                  gsl_interp_accel *accel,
                                                  bool angle = false)
    {
        optional<double> result;

        try
        {
            vector<double> y = util::optionals_to_values(yo);
            result = interpolate_vector_(x0, x, {y.data(), y.size()}, interp, accel);
        }
        catch (std::bad_optional_access)
        {
            result = std::nullopt;
        }

        return result;
    }

    Ephemeris::Datum Ephemeris::interpolate(const double mjd0) const
    {
        if ((mjd0 < data_.front().mjd) || (mjd0 > data_.back().mjd))
            throw std::runtime_error("Interpolation beyond ephemeris time range not supported");

        // order of polynomial interpolation is n - 1
        const gsl_interp_type *interp_type;
        int n = std::min(4, (int)data_.size());
        if (n < 2)
            throw std::runtime_error("Cannot interpolate with an ephemeris of length " + std::to_string(n));
        else if (n == 2)
            interp_type = gsl_interp_linear;
        else
            interp_type = gsl_interp_polynomial;

        // find the nearest segment; left and right are inclusive
        auto right = std::lower_bound(data_.cbegin(), data_.cend(), mjd0,
                                      [](const Datum &d, const double mjd)
                                      { return d.mjd.value() < mjd; });
        auto left = right - n + 1;

        // shift to center on mjd, taking care with integer math
        std::advance(right, (n - 2) / 2);
        std::advance(left, (n - 2) / 2);

        if (left < data_.cbegin())
        {
            int offset = std::distance(left, data_.cbegin());
            std::advance(right, offset);
            std::advance(left, offset);
        }

        if (right >= data_.cend())
        {
            int offset = std::distance(right, data_.cend() - 1);
            std::advance(right, offset);
            std::advance(left, offset);
        }

        if ((left < data_.cbegin()) || right >= data_.cend())
            throw std::runtime_error("Interpolation range is out of bounds.");

        // use the limits to define the segment that we will interpolate
        const Ephemeris segment(target_, {left, right + 1});
        const vector<double> mjd = util::optionals_to_values(segment.mjd());

        unique_interp_accel_ptr accel(gsl_interp_accel_alloc(), &gsl_interp_accel_free);
        unique_interp_ptr interp(gsl_interp_alloc(interp_type, n), &gsl_interp_free);

        Datum d;
        d.mjd = mjd0;
        d.tmtp = interpolate_optional_vector_(mjd0, mjd, segment.tmtp(), interp.get(), accel.get());
        d.mu = interpolate_optional_vector_(mjd0, mjd, segment.mu(), interp.get(), accel.get());
        d.mu_theta = interpolate_optional_vector_(mjd0, mjd, segment.mu_theta(), interp.get(), accel.get(), true);
        d.unc_a = interpolate_optional_vector_(mjd0, mjd, segment.unc_a(), interp.get(), accel.get());
        d.unc_b = interpolate_optional_vector_(mjd0, mjd, segment.unc_b(), interp.get(), accel.get());
        d.unc_theta = interpolate_optional_vector_(mjd0, mjd, segment.unc_theta(), interp.get(), accel.get(), true);
        d.rh = interpolate_optional_vector_(mjd0, mjd, segment.rh(), interp.get(), accel.get());
        d.delta = interpolate_optional_vector_(mjd0, mjd, segment.delta(), interp.get(), accel.get());
        d.phase = interpolate_optional_vector_(mjd0, mjd, segment.phase(), interp.get(), accel.get(), true);
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
            auto point = p.as_s2point();
            vx.push_back(point.x());
            vy.push_back(point.y());
            vz.push_back(point.z());
        }

        double x = interpolate_vector_(mjd0, mjd, {vx.data(), vx.size()}, interp.get(), accel.get());
        double y = interpolate_vector_(mjd0, mjd, {vy.data(), vy.size()}, interp.get(), accel.get());
        double z = interpolate_vector_(mjd0, mjd, {vz.data(), vz.size()}, interp.get(), accel.get());

        S2LatLng point(S2Point(x, y, z).Normalize());

        d.ra = point.lng().degrees();
        d.dec = point.lat().degrees();

        return d;
    }
}