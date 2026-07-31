#include "config.h"

#include "ephemeris.h"
#include "util/math.h"
#include "util/optional.h"

using std::cerr;
using std::endl;

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

    Ephemeris::Datum Ephemeris::interpolate(const double mjd) const
    {
        if ((mjd < data_.front().mjd) || (mjd > data_.back().mjd))
            throw std::runtime_error("Interpolation beyond ephemeris time range not supported");

        // find the nearest segment
        auto right = std::lower_bound(data_.cbegin(), data_.cend(), mjd,
                                      [](const Datum &d, const double mjd)
                                      { return d.mjd.value() < mjd; });
        auto left = right;
        if (right == data_.cbegin())
            right = std::next(right);
        else
            left = std::prev(left);

        const gsl_interp_type *interp_type;

        // order of polynomial is n - 1, we allow up to third order
        int n = std::min(4, (int)data_.size());
        if (n < 2)
            throw std::runtime_error("Cannot interpolate with an ephemeris of length " + std::to_string(n));

        // 2 point interpolation --> linear
        if (n == 2)
            interp_type = gsl_interp_linear;

        // 3 or 4 points --> polynomial
        if (n > 2)
        {
            interp_type = gsl_interp_polynomial;

            // expand to include next nearest segments, to the left is preferred
            if (left == data_.cbegin())
                right = std::next(right);
            else
                left = std::prev(left);
        }

        // 4 points
        if (n == 4)
        {
            // expand, to the right is preferred
            if (right == data_.cend() - 1)
                left = std::prev(left);
            else
                right = std::next(right);
        }

        // use the limits to define the segment that we will interpolate
        const Ephemeris segment(target_, {left, right + 1});
        const vector<double> mjd0 = util::optionals_to_values(segment.mjd());

        unique_interp_accel_ptr accel(gsl_interp_accel_alloc(), &gsl_interp_accel_free);
        unique_interp_ptr interp(gsl_interp_alloc(interp_type, n), &gsl_interp_free);

        Datum d;
        d.mjd = mjd;
        d.tmtp = interpolate_optional_vector_(mjd, mjd0, segment.tmtp(), interp.get(), accel.get());
        d.ra = interpolate_optional_vector_(mjd, mjd0, segment.ra(), interp.get(), accel.get());
        d.dec = interpolate_optional_vector_(mjd, mjd0, segment.dec(), interp.get(), accel.get());
        d.mu = interpolate_optional_vector_(mjd, mjd0, segment.mu(), interp.get(), accel.get());
        d.mu_theta = interpolate_optional_vector_(mjd, mjd0, segment.mu_theta(), interp.get(), accel.get());
        d.unc_a = interpolate_optional_vector_(mjd, mjd0, segment.unc_a(), interp.get(), accel.get());
        d.unc_b = interpolate_optional_vector_(mjd, mjd0, segment.unc_b(), interp.get(), accel.get());
        d.unc_theta = interpolate_optional_vector_(mjd, mjd0, segment.unc_theta(), interp.get(), accel.get());
        d.rh = interpolate_optional_vector_(mjd, mjd0, segment.rh(), interp.get(), accel.get());
        d.delta = interpolate_optional_vector_(mjd, mjd0, segment.delta(), interp.get(), accel.get());
        d.phase = interpolate_optional_vector_(mjd, mjd0, segment.phase(), interp.get(), accel.get());
        d.selong = interpolate_optional_vector_(mjd, mjd0, segment.selong(), interp.get(), accel.get());
        d.true_anomaly = interpolate_optional_vector_(mjd, mjd0, segment.true_anomaly(), interp.get(), accel.get());
        d.sangle = interpolate_optional_vector_(mjd, mjd0, segment.sangle(), interp.get(), accel.get());
        d.vangle = interpolate_optional_vector_(mjd, mjd0, segment.vangle(), interp.get(), accel.get());
        d.vmag = interpolate_optional_vector_(mjd, mjd0, segment.vmag(), interp.get(), accel.get());

        return d;
    }
}