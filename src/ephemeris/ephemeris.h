#ifndef SBS_EPHEMERIS_H_
#define SBS_EPHEMERIS_H_

#include <iostream>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include <boost/json.hpp>
#include <s2/s2point.h>
#include <s2/s2region_term_indexer.h>
#include <s2/s2region_union.h>

#include "config.h"
#include "date.h"
#include "moving_target.h"
#include "observatory.h"

using std::optional;
using std::string;
using std::tuple;
using std::unique_ptr;
using std::vector;
namespace json = boost::json;

namespace sbsearch::ephemeris
{
    class Ephemeris
    {
    public:
        // One entry in the ephemeris table.
        //
        // - mjd: modified Julian date, UTC
        // - tmtp: T-Tp, time from nearest perihelion (osculating elements), days
        // - sky coordinates, International Celestial Reference Frame:
        //   + ra, dec: in degrees
        //   + vertices: an S2Point vector (initialized via from_vertices)
        // - mu: proper motion, arcsec/minute
        // - mu_theta: position angle of proper motion, deg
        // - unc_a, unc_b: the semi-major and -minor axes of the uncertainty, arcsec
        // - unc_theta: the position angle of the uncertainty ellipse, deg
        //   semi-major axis, deg east oo north.
        // - rh: heliocentric distance, au
        // - delta: observer-target distance, au
        // - phase: sun-target-observer angle, deg
        // - selong: solar elongation, deg
        // - true_anomaly: just that, deg
        // - sangle: the position angle of the comet-sun vector, deg
        // - vangle: the position angle of the comet velocity vector, deg
        // - vmag: an approximate visual magnitude
        //
        // mjd, ra, dec, mu, mu_theta, rh, delta, phase are required.
        struct Datum
        {
            // time
            double mjd;
            optional<double> tmtp;
            // sky coordinates
            double ra;
            double dec;
            double mu;
            double mu_theta;
            optional<double> unc_a;
            optional<double> unc_b;
            optional<double> unc_theta;
            // geometry
            double rh;
            double delta;
            double phase;
            optional<double> selong;
            optional<double> true_anomaly;
            optional<double> sangle;
            optional<double> vangle;
            // other
            optional<double> vmag;

            // Get mjd as a date.
            inline Date date() const { return Date(mjd); };

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

            bool operator==(const Datum &other) const;
            bool operator!=(const Datum &other) const;
        };

        typedef vector<Datum> Data;

        // Output formatting.
        struct Format
        {
            Date::Format date = Date::Format::MJD;
        } format;

        // Initialize
        Ephemeris(const MovingTarget &target, const Data &data, const Format f = {Date::Format::MJD});

        // default constructor makes an empty ephemeris
        Ephemeris() : Ephemeris(MovingTarget(), {}) {};

        // copy constructor
        Ephemeris(const Ephemeris &other) = default;

        // move constructor
        Ephemeris(Ephemeris &&other) = default;

        // copy assignment
        Ephemeris &operator=(const Ephemeris &other) = default;

        // move assignment
        Ephemeris &operator=(Ephemeris &&other) = default;

        // validate ephemeris data
        bool isValid() const;

        // If the ephemeris is a single point, then values will be directly
        // printed with the ostream, otherwise the ephemeris will be printed as
        // a table.
        friend std::ostream &operator<<(std::ostream &os, const Ephemeris &ephemeris);

        // Return a single epoch (vertex) from the ephemeris.
        //
        // If `k<0`, then the index is relative to the end.
        Ephemeris operator[](const int k) const;

        // Return a slice of the ephemeris, from `start` up to `stop`, indexed by vertex.
        Ephemeris slice(const int start) const;
        Ephemeris slice(const int start, const int stop) const;

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
        vector<string> date() const; // "null" if mjd is null
        vector<optional<double>> tmtp() const;
        vector<double> ra() const;
        vector<double> dec() const;
        vector<double> mu() const;
        vector<double> mu_theta() const;
        vector<optional<double>> unc_a() const;
        vector<optional<double>> unc_b() const;
        vector<optional<double>> unc_theta() const;
        vector<double> rh() const;
        vector<double> delta() const;
        vector<double> phase() const;
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

        // Get ephemeris segment as an ephemeris object, if `k<0`, then the
        // index is relative to the end.
        Ephemeris segment(const int k) const;

        // Vector of ephemeris segments
        vector<Ephemeris> segments() const;

        // Append data.
        // mjd must follow in time.
        void append(const Datum &new_datum);
        void append(Datum &&new_datum);
        void append(const Data &new_data);
        void append(Data &&new_data);

        // Append the ephemeris.
        // Must have the same target and mjd must follow in time.
        void append(const Ephemeris &new_eph);
        void append(Ephemeris &&new_eph);

    private:
        int num_vertices_, num_segments_;
        MovingTarget target_;
        Data data_;

        // Get a vector of data member values.
        template <typename R, typename T>
        R get_data_vector(T Datum::*member) const;

        // index normalization maps the range (-max, max) to [0, max).
        static int normalize_index(const int k, const int max);

        // throws runtime_error if a.target() != b.target()
        static void check_target_id(const Ephemeris &a, const Ephemeris &b);

        // throws runtime_error if b does not follow a in time, and if b is not monotonically increasing
        static void check_relative_time(const Ephemeris &a, const Ephemeris::Data &b);
    };

    /** Implementation ***********************************************/

    template <typename R, typename T>
    R Ephemeris::get_data_vector(T Ephemeris::Datum::*member) const
    {
        R v(num_vertices_);
        std::transform(data_.begin(), data_.end(), v.begin(), std::mem_fn(member));
        return v;
    };
}

#endif // SBS_EPHEMERIS_H_
