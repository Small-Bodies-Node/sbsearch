#ifndef SBSEARCH_UTIL_SPHERICAL_H_
#define SBSEARCH_UTIL_SPHERICAL_H_

#include <vector>
#include <set>
#include <string>
#include <s2/s2latlng.h>
#include <s2/s2latlng_rect.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>

using std::string;
using std::vector;

namespace sbsearch::util
{
    // calculate the position angle (angle from North) from point a to point b
    double position_angle(const S2Point &a, const S2Point &b);

    // offset `distance` from `point` along `position_angle`.
    S2LatLng offset_by(const S2LatLng &point, const S1Angle &position_angle, const S1Angle &distance);

    // Radial distance from center of ellipse to point at angle phi.
    double ellipse_rho(const double a, const double b, const double phi);

    // Generate an approximate ellipse on the sphere, composed of n points, angles in radians.
    vector<S2LatLng> ellipse(const int n, const S2LatLng &center, const double &a, const double &b, const double &theta);

    // Points on the ellipse with tangent slopes aligned with position_angle.
    std::array<S2LatLng, 2> ellipse_tangent_points(const S2LatLng &center,
                                                   const double a,
                                                   const double b,
                                                   const double theta,
                                                   const double position_angle);
}

#endif // SBSEARCH_UTIL_SPHERICAL_H_