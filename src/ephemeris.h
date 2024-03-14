#ifndef SBS_EPHEMERIS_H_
#define SBS_EPHEMERIS_H_

#include "config.h"

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <s2/s2point.h>
#include <s2/s2polyline.h>
#include <s2/s2region_term_indexer.h>

#include "moving_target.h"
#include "observatory.h"

#define UNDEF_TIME -1
#define UNDEF_ANGLE -999
#define UNDEF_UNC -1

using std::string;
using std::tuple;
using std::unique_ptr;
using std::vector;

namespace sbsearch
{
    class Ephemeris
    {
    public:
        // One entry in the ephemeris table.
        // - mjd: modified Julian date, UTC
        // - tmtp: T-Tp, time from nearest perihelion (osculating elements), days
        // - sky coordinates, International Celestial Reference Frame:
        //   + ra, dec: in degrees
        //   + vertices: an S2Point vector (initialized via from_vertices)
        // - unc_a, unc_b: the semi-major and -minor axes of the uncertainty
        // - unc_theta: the position angle of the uncertainty ellipse
        //   semi-major axis, deg east of north.
        // - rh: heliocentric distance, au
        // - delta: observer-target distance, au
        // - phase: sun-target-observer angle, deg
        // - selong: solar elongation, deg
        // - true_anomaly: just that, deg
        // - sangle: the position angle of the comet-sun vector, deg
        // - vangle: the position angle of the comet velocity vector, deg
        // - vmag: an approximate visual magnitude
        struct Datum
        {
            // time
            double mjd = UNDEF_TIME;
            double tmtp = UNDEF_TIME;
            // sky coordinates
            double ra = UNDEF_ANGLE;
            double dec = UNDEF_ANGLE;
            double unc_a = UNDEF_UNC;
            double unc_b = UNDEF_UNC;
            double unc_theta = 0;
            // geometry
            double rh = 0;
            double delta = 0;
            double phase = UNDEF_ANGLE;
            double selong = UNDEF_ANGLE;
            double true_anomaly = UNDEF_ANGLE;
            double sangle = UNDEF_ANGLE;
            double vangle = UNDEF_ANGLE;
            // other
            double vmag = 99;

            // set RA, Dec from S2LatLng
            void radec(const S2LatLng &ll)
            {
                ra = ll.lng().degrees();
                dec = ll.lat().degrees();
            }

            // set RA, Dec from S2Point
            void radec(const S2Point &point)
            {
                radec(S2LatLng(point));
            }

            // RA, Dec as S2LatLng
            S2LatLng as_s2latlng() const
            {
                return S2LatLng::FromDegrees(dec, ra).Normalized();
            }

            // RA, Dec as S2Point
            S2Point as_s2point() const
            {
                return as_s2latlng().ToPoint();
            }

            bool operator==(const Datum &other) const;
            bool operator!=(const Datum &other) const;
        };

        typedef vector<Datum> Data;

        // For ephemeris extrapolation: BACKWARDS to extrapolate before the
        // first vertex, FORWARDS to extrapolate beyond the last vertex.
        enum class Extrapolate : uint8
        {
            BACKWARDS,
            FORWARDS
        };

        // Initialize
        Ephemeris(const MovingTarget target, Data data);

        // default constructor makes an empty ephemeris
        Ephemeris() : Ephemeris(MovingTarget(), {}){};

        // validate ephemeris data
        bool isValid() const;

        // Ephemeris search options: may use uncertainties, padding, or both.
        struct Options
        {
            bool use_uncertainty = false;
            double padding = 0; // arcsec

            // true if padding or use_uncertainty are enabled
            bool padding_enabled() const
            {
                return use_uncertainty | (padding > 0);
            }
        };

        // options, may be changed at any time
        inline const Options &options() const { return options_; }
        inline Options *mutable_options() { return &options_; }

        // output
        //
        // Format options; zero for default.
        struct Format
        {
            size_t designation_width = 0;
            size_t moving_target_id_width = 0;
            size_t tmtp_width = 0;
        } format;

        // Calculate column widths for stream output.
        Format format_widths() const;

        // If the ephemeris is a single point, then it will be printed without a
        // terminating new-line, otherwise the ephemeris will be printed as a table.
        friend std::ostream &operator<<(std::ostream &os, const Ephemeris &ephemeris);

        // Return a single epoch from the ephemeris.
        //
        // If `k<0`, then the index is relative to the end.
        const Ephemeris operator[](const int k) const;

        // Return a slice of the ephemeris.
        const Ephemeris slice(const int start);
        const Ephemeris slice(const int start, const int stop);

        // equality tests
        bool operator==(const Ephemeris &other) const;
        bool operator!=(const Ephemeris &other) const { return !((*this) == other); };
        // Number of ephemeris vertices
        int num_vertices() const;

        // Property getters
        inline const MovingTarget &target() const { return target_; }
        inline const Data &data() const { return data_; };
        // If `k<0`, then the index is relative to the end.
        const Datum &data(const int k) const;

        // Array access
        vector<double> mjd() const;
        vector<double> tmtp() const;
        vector<double> ra() const;
        vector<double> dec() const;
        vector<double> unc_a() const;
        vector<double> unc_b() const;
        vector<double> unc_theta() const;
        vector<double> rh() const;
        vector<double> delta() const;
        vector<double> phase() const;
        vector<double> selong() const;
        vector<double> true_anomaly() const;
        vector<double> sangle() const;
        vector<double> vangle() const;
        vector<double> vmag() const;

        // RA, Dec as S2Point(s)
        S2Point vertex(const int k) const;
        vector<S2Point> vertices() const;

        // target is mutable
        void target(const MovingTarget &new_target) { target_ = new_target; }

        // Number of ephemeris segments
        int num_segments() const;

        // Append the data.
        // mjd must follow in time.
        void append(const Data &new_data);

        // Append the ephemeris.
        // Must have the same target and mjd must follow in time.
        void append(const Ephemeris &eph);

        // Get ephemeris segment as an ephemeris object, if `k<0`, then the
        // index is relative to the end.
        Ephemeris segment(const int k) const;

        // Vector of ephemeris segments
        vector<Ephemeris> segments() const;

        // Ephemeris as a polyline
        S2Polyline as_polyline() const;

        // Offset the ephemeris for parallax.
        Ephemeris parallax_offset(const Observatory &observatory);

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
        void pad(const vector<double> &para, const vector<double> &perp, S2Polygon &polygon) const;
        void pad(const double para, const double perp, S2Polygon &polygon) const;

        // Pad a region around the ephemeris.
        //
        // The ephemeris is extended by vector `a` along direction `theta`, and
        // by vector `b` along `theta + 90 deg`.  Essentially a quadrilateral
        // approximation to a series of ellipses.
        //
        /// `a` and `b` in units of arcsec, `theta` in units of degrees east of north
        void pad(const double a, const double b, const double theta, S2Polygon &polygon) const;
        void pad(const vector<double> &a, const vector<double> &b, const vector<double> &theta, S2Polygon &polygon) const;

        //
        void as_polygon(S2Polygon &polygon) const;

    private:
        int num_vertices_, num_segments_;
        MovingTarget target_;
        Data data_;
        Options options_;
        int normalize_index(const int i, const int max) const;
    };
}

#endif // SBS_EPHEMERIS_H_