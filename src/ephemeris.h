#ifndef SBS_EPHEMERIS_H_
#define SBS_EPHEMERIS_H_

#include "config.h"

#include <string>
#include <vector>
#include <memory>
#include <s2/s2point.h>
#include <s2/s2polyline.h>
#include <s2/s2region_term_indexer.h>

#define UNDEF_UNC -1
#define UNDEF_OBJECT_ID -1

using std::string;
using std::unique_ptr;
using std::vector;

namespace sbsearch
{
    class Ephemeris
    {
    public:
        // For ephemeris extrapolation: BACKWARDS to extrapolate before the
        // first vertex, FORWARDS to extrapolate beyond the last vertex.
        enum class Extrapolate : uint8
        {
            BACKWARDS,
            FORWARDS
        };

        // Ephemeris search options: may use uncertainties, padding, or both.
        struct Options
        {
            bool use_uncertainty = false;
            double padding = 0; // radians

            // true if padding or use_uncertainty are enabled
            bool padding_enabled() const
            {
                return use_uncertainty | (padding > 0);
            }
        };

        // Initialize from vectors
        // - object_id is a unique object identifier for this ephemeris
        // - vertices are RA and Dec, International Celestial Reference Frame
        // - mjd is modified Julian date, UTC
        // - rh is heliocentric distance, au
        // - delta is observer-target distance, au
        // - phase is sun-target-observer angle, deg
        // - unc_a, unc_b are the semi-major and -minor axes of the uncertainty
        //   ellipse
        // - unc_theta is the position angle of the uncertainty ellipse
        //   semi-major axis, deg east of north.
        Ephemeris(const int object_id,
                  const vector<S2Point> &vertices,
                  const vector<double> &mjd,
                  const vector<double> &rh,
                  const vector<double> &delta,
                  const vector<double> &phase,
                  const vector<double> &unc_a,
                  const vector<double> &unc_b,
                  const vector<double> &unc_theta);

        // Convenience function, mostly for testing
        Ephemeris(const int object_id,
                  const vector<S2Point> &vertices,
                  const vector<double> &mjd,
                  const vector<double> &rh,
                  const vector<double> &delta,
                  const vector<double> &phase)
            : Ephemeris(object_id, vertices, mjd, rh, delta, phase,
                        vector<double>(vertices.size(), UNDEF_UNC),
                        vector<double>(vertices.size(), UNDEF_UNC),
                        vector<double>(vertices.size(), UNDEF_UNC)) {}

        // object_id is mutable
        void object_id(int new_object_id) { object_id_ = new_object_id; }

        // default constructor makes an empty ephemeris
        Ephemeris() : Ephemeris(UNDEF_OBJECT_ID, {}, {}, {}, {}, {}, {}, {}, {}){};

        // return a single-point ephemeris, if `k<0`, then the index is relative
        // to the end.
        Ephemeris operator[](const int k) const;

        // validate ephemeris data
        bool isValid() const;

        // output
        //
        // Format options; zero for default.
        struct Format
        {
            size_t object_id_width = 0;
        } format;

        // If the ephemeris is a single point, then it will be printed without a
        // terminating new-line, otherwise the ephemeris will be printed as a table.
        friend std::ostream &operator<<(std::ostream &os, const Ephemeris &ephemeris);

        // equality tests
        bool is_equal(const Ephemeris &other) const;

        // options, may be changed at any time
        inline const Options &options() const { return options_; }
        inline Options *mutable_options() { return &options_; }

        // Number of ephemeris vertices
        int num_vertices() const;

        // Property getters, if `k<0`, then the index is relative to the end.
        inline const int &object_id() const { return object_id_; }
        inline const vector<S2Point> &vertices() const { return vertices_; }
        inline const vector<double> &mjd() const { return mjd_; }
        inline const vector<double> &rh() const { return rh_; }
        inline const vector<double> &delta() const { return delta_; }
        inline const vector<double> &phase() const { return phase_; }
        inline const vector<double> &unc_a() const { return unc_a_; }
        inline const vector<double> &unc_b() const { return unc_b_; }
        inline const vector<double> &unc_theta() const { return unc_theta_; }

        const S2Point &vertex(const int k) const;
        inline const double &mjd(const int k) const { return getter(mjd_, k); }
        inline const double &rh(const int k) const { return getter(rh_, k); }
        inline const double &delta(const int k) const { return getter(delta_, k); }
        inline const double &phase(const int k) const { return getter(phase_, k); }
        inline const double &unc_a(const int k) const { return getter(unc_a_, k); }
        inline const double &unc_b(const int k) const { return getter(unc_b_, k); }
        inline const double &unc_theta(const int k) const { return getter(unc_theta_, k); }

        // vertex as RA, Dec
        double ra(const int k) const;
        double dec(const int k) const;

        // Number of ephemeris segments
        int num_segments() const;

        // Append the ephemeris, must have the same object_id.
        void append(const Ephemeris &eph);

        // Get ephemeris segment as an ephemeris object, if `k<0`, then the
        // index is relative to the end.
        Ephemeris segment(const int k) const;

        // Vector of ephemeris segments
        vector<Ephemeris> segments() const;

        // Ephemeris as a polyline
        S2Polyline as_polyline() const;

        // Linearly interpolate ephemeris to time `mjd0`.
        Ephemeris interpolate(const double mjd0) const;

        // Linearly (on the sphere) extrapolate ephemeris by amount `distance`
        // in radians
        Ephemeris extrapolate(const double distance, Extrapolate direction) const;

        /* Get a subsample of the ephemeris based on the given date range

        Comet and asteroid motion is non-linear, but this method uses a linear
        approximation.  (Non-linearity should be addressed with finer ephemeris
        steps.)  In order to minimize errors, only test the nearest segment(s)
        to the observation.  For example:

              0               1
        |----------|--------------------|

        Segment 0: t0 = 0 dt = 1 da = 10 deg

        Segment 1: t0 = 1 dt = 1 da = 20 deg

        Average proper motion: 30 deg / 2 days = 15 deg / day

        Linear interpolation to t = 1? --> 15 deg

        But we wanted 10 deg.

        */
        Ephemeris subsample(const double mjd_start, const double mjd_stop) const;

        // Pad a region around the ephemeris and return the result as a polygon.
        // `para` is padding parallel to the ephemeris, `perp` is perpendicular
        // to it.  Both values in radians.  For vectors, there must be one
        // element per ephemeris vertex.  The offsets must be less than 90 deg.
        // `para` and `perp` in units of arcsec.
        S2Polygon pad(const vector<double> &para, const vector<double> &perp) const;
        S2Polygon pad(const double para, const double perp) const;

        // Pad a region around the ephemeris.
        //
        // The ephemeris is extended by vector `a` along direction `theta`, and
        // by vector `b` along `theta + 90 deg`.  Essentially a quadrilateral
        // approximation to a series of ellipses.
        //
        /// `a` and `b` in units of arcsec, `theta` in units of degrees east of north
        S2Polygon pad(const double a, const double b, const double theta) const;
        S2Polygon pad(const vector<double> &a, const vector<double> &b, const vector<double> &theta) const;

        //
        S2Polygon as_polygon() const;

    private:
        int num_vertices_, num_segments_, object_id_;
        vector<S2Point> vertices_;
        vector<double> mjd_, rh_, delta_, phase_, unc_a_, unc_b_, unc_theta_;
        Options options_;

        // get single value
        const double &getter(const vector<double> &property, const int k) const;
    };
}

#endif // SBS_EPHEMERIS_H_