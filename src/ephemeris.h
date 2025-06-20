#ifndef SBS_EPHEMERIS_H_
#define SBS_EPHEMERIS_H_

#include "config.h"

#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include <boost/json.hpp>
#include <s2/s2point.h>
#include <s2/s2polyline.h>
#include <s2/s2region_term_indexer.h>
#include <s2/s2region_union.h>

#include "moving_target.h"
#include "observatory.h"

using std::optional;
using std::string;
using std::tuple;
using std::unique_ptr;
using std::vector;
namespace json = boost::json;

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
        // - mu: proper motion, arcsec/minute
        // - mu_theta: position angle of proper motion, deg
        // - unc_a, unc_b: the semi-major and -minor axes of the uncertainty, arcsec
        // - unc_theta: the position angle of the uncertainty ellipse, deg
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
            optional<double> mjd;
            optional<double> tmtp;
            // sky coordinates
            optional<double> ra;
            optional<double> dec;
            optional<double> mu;
            optional<double> mu_theta;
            optional<double> unc_a;
            optional<double> unc_b;
            optional<double> unc_theta;
            // geometry
            optional<double> rh;
            optional<double> delta;
            optional<double> phase;
            optional<double> selong;
            optional<double> true_anomaly;
            optional<double> sangle;
            optional<double> vangle;
            // other
            optional<double> vmag;

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
                return S2LatLng::FromDegrees(dec.value(), ra.value()).Normalized();
            }

            // RA, Dec as S2Point
            S2Point as_s2point() const
            {
                return as_s2latlng().ToPoint();
            }

            bool operator==(const Datum &other) const;
            bool operator!=(const Datum &other) const;

            // Return data as JSON object
            json::object as_json();
        };

        typedef vector<Datum> Data;

        // For ephemeris extrapolation: BACKWARDS to extrapolate before the
        // first vertex, FORWARDS to extrapolate beyond the last vertex.
        enum struct Extrapolate : uint8
        {
            BACKWARDS,
            FORWARDS
        };

        // Initialize
        Ephemeris(const MovingTarget target, Data data);

        // default constructor makes an empty ephemeris
        Ephemeris() : Ephemeris(MovingTarget(), {}) {};

        // validate ephemeris data
        bool isValid() const;

        // Ephemeris options
        struct Options
        {
            bool use_uncertainty = false; // increase search area using uncertainty
        };

        // options, may be changed at any time
        inline const Options &options() const { return options_; }
        inline Options *mutable_options() { return &options_; }

        // Stream output format options.
        struct Format
        {
            // Display dates as "mjd" or "calendar"?
            enum struct DateFormat : uint8
            {
                MJD,
                CALENDAR
            };
            DateFormat date = DateFormat::MJD;
        } format;

        // If the ephemeris is a single point, then values will be directly
        // printed with the ostream, otherwise the ephemeris will be printed as
        // a table.
        friend std::ostream &operator<<(std::ostream &os, const Ephemeris &ephemeris);

        // Return a single epoch (vertex) from the ephemeris.
        //
        // If `k<0`, then the index is relative to the end.
        const Ephemeris operator[](const int k) const;

        // Return a slice of the ephemeris, from `start` up to `stop`, indexed by vertex.
        const Ephemeris slice(const int start) const;
        const Ephemeris slice(const int start, const int stop) const;

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
        vector<optional<double>> mjd() const;
        vector<string> date() const; // "null" if mjd is null
        vector<optional<double>> tmtp() const;
        vector<optional<double>> ra() const;
        vector<optional<double>> dec() const;
        vector<optional<double>> mu() const;
        vector<optional<double>> mu_theta() const;
        vector<optional<double>> unc_a() const;
        vector<optional<double>> unc_b() const;
        vector<optional<double>> unc_theta() const;
        vector<optional<double>> rh() const;
        vector<optional<double>> delta() const;
        vector<optional<double>> phase() const;
        vector<optional<double>> selong() const;
        vector<optional<double>> true_anomaly() const;
        vector<optional<double>> sangle() const;
        vector<optional<double>> vangle() const;
        vector<optional<double>> vmag() const;

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

        // Split ephemeris in segments of approximate length `length` in degrees and `time` in days.
        vector<Ephemeris> split(double length, double time) const;

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
        steps.)  In order to minimize errors, only test the segment(s) nearest
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

        // Convert the ephemeris into a vector of polygons, with optional
        // padding in arcmin.
        //
        // The area will depend on the `use_uncertainty` option, but the padding
        // around the ephemeris will be at least 2", based on results in the testing suite.
        //
        // The ephemeris is described by connecting parallelograms that
        // circumscribe ellipses with semi-major axes `a`, semi-minor axes `b`,
        // with `a` aligned along angle `theta` (E of N).
        vector<unique_ptr<S2Polygon>> as_polygons(double padding = 0) const;

        // Return data as JSON array
        json::array as_json();

    private:
        int num_vertices_, num_segments_;
        MovingTarget target_;
        Data data_;
        Options options_;
        int normalize_index(const int i, const int max) const;
    };

    std::istream &operator>>(std::istream &in, Ephemeris::Format::DateFormat &date_format);
}

#endif // SBS_EPHEMERIS_H_