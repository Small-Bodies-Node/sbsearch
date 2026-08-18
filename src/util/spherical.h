#ifndef SBSEARCH_UTIL_SPHERICAL_H_
#define SBSEARCH_UTIL_SPHERICAL_H_

#include <memory>
#include <set>
#include <string>
#include <vector>
#include <s2/s2latlng.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>

using std::string;
using std::vector;

namespace sbsearch::util
{
    // Archaversine
    double ahav(const double x);

    // Circumscribe ellipse with a parallelogram aligned with `position_angle`.
    // All angles in radians.
    std::unique_ptr<S2Polygon> circumscribe_ellipse(const S2LatLng &center,
                                                    const double a,
                                                    const double b,
                                                    const double theta,
                                                    const double position_angle);

    // Generate an approximate ellipse on the sphere, composed of n points,
    // angles in radians.
    vector<S2LatLng> ellipse(const int n,
                             const S2LatLng &center,
                             const double &a,
                             const double &b,
                             const double &theta);

    // Radial distance from center of ellipse to point at angle phi, all angles
    // in radians.
    double ellipse_rho(const double a, const double b, const double phi);

    // Half-width of an ellipse along `position_angle`.  For example, the
    // half-width is a when position_angle == a, b when position_angle == 90
    // deg.
    double ellipse_half_width(const double a,
                              const double b,
                              const double theta,
                              const double position_angle);

    // Points on the ellipse with tangent slopes aligned with position_angle,
    // angles in radians.  The point CCW from theta will be first.
    std::array<S2LatLng, 2> ellipse_tangent_points(const S2LatLng &center,
                                                   const double a,
                                                   const double b,
                                                   const double theta,
                                                   const double position_angle);

    // Haversine
    double hav(const double theta);

    // The smallest angle on the circle between angles a and b.
    // S1Angle inner_angle(const S1Angle &a, const S1Angle &b);

    // offset `distance` from `point` along `position_angle`.
    S2LatLng offset_by(const S2LatLng &point, const S1Angle &position_angle, const S1Angle &distance);

    // calculate the position angle (angle from North) from point a to point b
    double position_angle(const S2Point &a, const S2Point &b);
    double position_angle(const S2LatLng &a, const S2LatLng &b);
}

#endif // SBSEARCH_UTIL_SPHERICAL_H_