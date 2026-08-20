#ifndef SBS_EPHEMERIS_INTERPOLATE_H_
#define SBS_EPHEMERIS_INTERPOLATE_H_

#include <optional>
#include <valarray>
#include <vector>
#include <gsl/gsl_interp.h>

#include "ephemeris/ephemeris.h"

using sbsearch::ephemeris::Ephemeris;
using std::optional;
using std::valarray;
using std::vector;

namespace sbsearch::ephemeris
{
    // Interpolate ephemeris to time `mjd0`.
    Ephemeris::Datum interpolate(const Ephemeris::Data &data, const double mjd0);
    Ephemeris interpolate(const Ephemeris &eph, const double mjd0);

    // Interpolate a vector as a function of x to a specific date.
    double interpolate_valarray_(const double x0,
                                 const vector<double> &x,
                                 const valarray<double> &y,
                                 gsl_interp *interp,
                                 gsl_interp_accel *accel,
                                 bool angle = false);

    // Interpolate a vector as a function of x to a specific date.
    double interpolate_vector_(const double x0,
                               const vector<double> &x,
                               const vector<double> &y,
                               gsl_interp *interp,
                               gsl_interp_accel *accel,
                               bool angle = false);

    // RAII with GNU Scientific Library pointers
    using unique_interp_accel_ptr = std::unique_ptr<gsl_interp_accel, decltype(&gsl_interp_accel_free)>;
    using unique_interp_ptr = std::unique_ptr<gsl_interp, decltype(&gsl_interp_free)>;

    // Interpolate a vector of optional values, and return nullopt if any values
    // are not defined.
    optional<double> interpolate_optional_vector_(const double x0,
                                                  const vector<double> &x,
                                                  const vector<optional<double>> &yo,
                                                  gsl_interp *interp,
                                                  gsl_interp_accel *accel,
                                                  bool angle = false);
}

#endif