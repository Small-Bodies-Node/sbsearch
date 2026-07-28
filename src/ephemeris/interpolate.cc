#include "config.h"

#include "ephemeris.h"
#include "util/math.h"
#include "util/optional.h"

using std::cerr;

// Support RAII with GNU Scientific Library pointers
using unique_interp_accel_ptr = std::unique_ptr<gsl_interp_accel, decltype(&gsl_interp_accel_free)>;
using unique_interp_ptr = std::unique_ptr<gsl_interp, decltype(&gsl_interp_free)>;

namespace sbsearch
{
    // Interpolation helper that handles vectors of optional values, and returns nullopt if any values are not defined.
    optional<double> interpolate_optional_vector_(const double mjd,
                                                  const vector<double> &mjd0,
                                                  const vector<optional<double>> &yo,
                                                  gsl_interp *interp,
                                                  gsl_interp_accel *accel)
    {
        optional<double> result;

        try
        {
            vector<double> y = util::optionals_to_values(yo);
            gsl_interp_init(interp, mjd0.data(), y.data(), y.size());
            result = gsl_interp_eval(interp, mjd0.data(), y.data(), mjd, accel);
        }
        catch (std::bad_optional_access)
        {
            result = std::nullopt;
        }

        return result;
    };

    Ephemeris Ephemeris::interpolate(const double mjd) const
    {
        if ((mjd < data_.front().mjd) || (mjd > data_.back().mjd))
            throw std::runtime_error("Interpolation beyond ephemeris time range not supported");

        // find the nearest segment
        auto right = std::lower_bound(data_.begin(), data_.end(), mjd,
                                      [](const Datum &d, const double mjd)
                                      { return d.mjd.value() < mjd; });
        auto left = std::prev(right);

        // return value
        Datum d;

        const gsl_interp_type *interp_type;

        // order of polynomial is n - 1, we allow up to third order
        int n = std::min(4, (int)data_.size());
        if (n == 2)
        {
            // 2 point interpolation --> linear
            interp_type = gsl_interp_linear;
        }

        if (n > 2)
        {
            // 3 or 4 points --> polynomial
            interp_type = gsl_interp_polynomial;

            // expand to include next nearest segments
            if (std::distance(data_.begin(), left) > 0)
                left = std::prev(left);
            else
                right = std::next(right);
        }

        if (n == 4)
        {
            // 4 points
            if (std::distance(right, data_.end()) > 1)
                right = std::next(right);
            else
                left = std::prev(left);
        }

        const Ephemeris eph(target_, {left, right + 1});
        const vector<double> mjd0 = util::optionals_to_values(eph.mjd());

        unique_interp_accel_ptr accel(gsl_interp_accel_alloc(), &gsl_interp_accel_free);
        unique_interp_ptr interp(gsl_interp_alloc(interp_type, n), &gsl_interp_free);

        d.mjd = mjd;
        d.tmtp = interpolate_optional_vector_(mjd, mjd0, eph.tmtp(), interp.get(), accel.get());
        d.ra = interpolate_optional_vector_(mjd, mjd0, eph.ra(), interp.get(), accel.get());
        d.dec = interpolate_optional_vector_(mjd, mjd0, eph.dec(), interp.get(), accel.get());
        d.mu = interpolate_optional_vector_(mjd, mjd0, eph.mu(), interp.get(), accel.get());
        d.mu_theta = interpolate_optional_vector_(mjd, mjd0, eph.mu_theta(), interp.get(), accel.get());
        d.unc_a = interpolate_optional_vector_(mjd, mjd0, eph.unc_a(), interp.get(), accel.get());
        d.unc_b = interpolate_optional_vector_(mjd, mjd0, eph.unc_b(), interp.get(), accel.get());
        d.unc_theta = interpolate_optional_vector_(mjd, mjd0, eph.unc_theta(), interp.get(), accel.get());
        d.rh = interpolate_optional_vector_(mjd, mjd0, eph.rh(), interp.get(), accel.get());
        d.delta = interpolate_optional_vector_(mjd, mjd0, eph.delta(), interp.get(), accel.get());
        d.phase = interpolate_optional_vector_(mjd, mjd0, eph.phase(), interp.get(), accel.get());
        d.selong = interpolate_optional_vector_(mjd, mjd0, eph.selong(), interp.get(), accel.get());
        d.true_anomaly = interpolate_optional_vector_(mjd, mjd0, eph.true_anomaly(), interp.get(), accel.get());
        d.sangle = interpolate_optional_vector_(mjd, mjd0, eph.sangle(), interp.get(), accel.get());
        d.vangle = interpolate_optional_vector_(mjd, mjd0, eph.vangle(), interp.get(), accel.get());
        d.vmag = interpolate_optional_vector_(mjd, mjd0, eph.vmag(), interp.get(), accel.get());

        return Ephemeris{target_, {d}};
    }
}