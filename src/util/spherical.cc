#include "config.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <ostream>
#include <memory>
#include <string>
#include <vector>
#include <s2/s1angle.h>
#include <s2/s2builder.h>
#include <s2/s2buffer_operation.h>
#include <s2/s2builder.h>
#include <s2/s2builderutil_lax_polygon_layer.h>
#include <s2/s2builderutil_s2polygon_layer.h>
#include <s2/s2debug.h>
#include <s2/s2error.h>
#include <s2/s2latlng.h>
#include <s2/s2lax_polygon_shape.h>
#include <s2/s2loop.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>
#include <sqlite3.h>

#include "spherical.h"
#include "../constants.h"
#include "../logging.h"
#include "../sofa/sofa.h"

using std::atan2;
using std::ceil;
using std::cos;
using std::floor;
using std::sin;
using std::string;
using std::tan;
using std::vector;

namespace sbsearch::util
{
    double position_angle(const S2Point &a, const S2Point &b)
    {
        // Meeus 1998, Astronomical Algorithms, p116
        S2LatLng aa(a);
        S2LatLng bb(b);
        S1Angle dra = bb.lng() - aa.lng();
        return atan2(sin(dra), cos(aa.lat()) * tan(bb.lat()) - sin(aa.lat()) * cos(dra));
    }

    S2LatLng offset_by(const S2LatLng &origin, const S1Angle &position_angle, const S1Angle &distance)
    {
        // Spherical trig based on astropy.coordinates.angle_utilities.offset_by()

        // Spherical trigonometry with:
        // Triangle:
        //   A: north pole (or really-really close to it)
        //   B: point
        //   C: result
        // With angles:
        //   A: change in longitude
        //   B: position angle
        //   C:
        // And sides:
        //   a: distance
        //   b: final co-latitude
        //   c: starting co-latitude

        double cos_a = cos(distance);
        double sin_a = sin(distance);
        double cos_c = sin(origin.lat());
        double sin_c = cos(origin.lat());
        double cos_B = cos(position_angle);
        double sin_B = sin(position_angle);

        double cos_b = cos_c * cos_a + sin_c * sin_a * cos_B;
        double xsin_A = sin_a * sin_B * sin_c;
        double xcos_A = cos_a - cos_b * cos_c;

        S1Angle A = S1Angle::Radians(std::atan2(xsin_A, xcos_A));
        bool small_sin_c = sin_c < 1e-12;
        if (small_sin_c)
        {
            A = S1Angle::Radians(PI_2 + cos_c * (PI_2 - position_angle.radians()));
        }

        S1Angle lon = origin.lng() + A;
        S1Angle lat = S1Angle::Radians(std::asin(cos_b));

        return S2LatLng(lat, lon);
    }

    double ellipse_rho(const double a, const double b, const double phi)
    {
        return a * b / std::sqrt(std::pow(b * cos(phi), 2) + std::pow(a * sin(phi), 2));
    }

    vector<S2LatLng> ellipse(const int n, const S2LatLng &center, const double &a, const double &b, const double &theta)
    {
        // all angles in radians
        vector<S2LatLng> e;
        const S1Angle th = S1Angle::Radians(theta);

        assert(n >= 4);
        for (int i = 0; i < n; i++)
        {
            const S1Angle phi = S1Angle::Radians(2 * PI * i / n);
            const S1Angle rho = S1Angle::Radians(ellipse_rho(a, b, phi.radians()));
            e.emplace_back(offset_by(center, th + phi, rho));
        }

        return e;
    }

    std::array<S2LatLng, 2> ellipse_tangent_points(const S2LatLng &center,
                                                   const double a,
                                                   const double b,
                                                   const double theta,
                                                   const double position_angle)
    {
        // the slope, transformed for an ellipse centered at 0, 0, with theta = 0
        const double m = std::tan(position_angle - theta);

        // tangent points are along pa and pa + 180
        const double u = std::sqrt(m * m * a * a + b * b);
        const double phi = std::atan2(-m * a * a / u, b * b / u);
        const double rho = ellipse_rho(a, b, phi);

        // offset center point to each tangent point
        const S2LatLng ll1 = offset_by(center, S1Angle::Radians(phi + theta), S1Angle::Radians(rho));
        const S2LatLng ll2 = offset_by(center, S1Angle::Radians(phi + theta + M_PI), S1Angle::Radians(rho));

        // the point to the left (CCW) will always be first
        S2Point mu = offset_by(center, S1Angle::Radians(position_angle), S1Angle::Degrees(1 * ARCSEC)).ToPoint() - center.ToPoint();
        S2Point p1 = ll1.ToPoint();
        double sinth1 = mu.x() * p1.y() + mu.y() * p1.z() + mu.z() * p1.x() - p1.y() * mu.z() - p1.z() * mu.x() - p1.x() * mu.y();

        if (sinth1 >= 0)
            return {ll1, ll2};
        else
            return {ll2, ll1};
    }
}