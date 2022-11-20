#ifndef EPHEMERIS_H_
#define EPHEMERIS_H_

#include <string>
#include <vector>
#include <memory>
#include <s2/s2point.h>
#include <s2/s2polyline.h>
#include <s2/s2region_term_indexer.h>

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

        // Initialize from vectors of vertices and times
        // - vertices are RA and Dec
        // - time is modified Julian date
        Ephemeris(const vector<S2Point> vertices, const vector<double> time);

        // validate ephemeris data
        bool isValid();

        // Number of ephemeris vertices
        int num_vertices();

        // Get vertex, if `k<0`, then the index is relative to the end.
        S2Point vertex(const int k);

        // Get vector of vertices
        vector<S2Point> vertices();

        // Get time, if `k<0`, then the index is relative to the end.
        double time(const int k);

        // Get vector of time
        vector<double> times();

        // Number of ephemeris segments
        int num_segments();

        // Get ephemeris segment as an ephemeris object, if `k<0`, then the
        // index is relative to the end.
        Ephemeris segment(const int k);

        // Vector of ephemeris segments
        vector<Ephemeris> segments();

        // Ephemeris as a polyline
        S2Polyline as_polyline();

        // Linearly (on the sphere) interpolate ephemeris to time `mjd`
        S2Point interpolate(const double mjd);

        // Linearly (on the sphere) extrapolate ephemeris by amount `distance`
        // in radians
        S2Point extrapolate(const double distance, Extrapolate direction);

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
        Ephemeris subsample(const double mjd_start, const double mjd_stop);

        // Pad a region around the ephemeris and return the result as a polygon.
        // `para` is padding parallel to the ephemeris, `perp` is perpendicular
        // to it.  Both values in radians.  For vectors, there must be one
        // element per ephemeris vertex.  The offsets must be less than 90 deg.
        unique_ptr<S2Polygon> pad(const vector<double> &para, const vector<double> &perp);
        unique_ptr<S2Polygon> pad(const double para, const double perp);

        // Generate ephemeris query terms
        vector<string> query_terms(S2RegionTermIndexer &indexer);

        // Query terms for a padded ephemeris
        vector<string> query_terms(S2RegionTermIndexer &indexer, const vector<double> &para, vector<double> &perp);
        vector<string> query_terms(S2RegionTermIndexer &indexer, const double para, const double perp);

    private:
        int num_vertices_, num_segments_;
        vector<S2Point> vertices_;
        vector<double> times_;
    };
}

#endif // EPHEMERIS_H_