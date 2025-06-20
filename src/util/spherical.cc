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

#include "polygon.h"
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
    double ahav(const double x)
    {
        return 2 * std::asin(std::sqrt(x));
    }

    std::unique_ptr<S2Polygon> circumscribe_ellipse(const S2LatLng &center,
                                                    const double a,
                                                    const double b,
                                                    const double theta,
                                                    const double position_angle)
    {
        auto [t1, t2] = ellipse_tangent_points(center, a, b, theta, position_angle);
        const S1Angle d = S1Angle::Radians(ellipse_rho(a, b, position_angle - theta));
        const S1Angle pa = S1Angle::Radians(position_angle);

        std::vector<S2Point> vertices({
            offset_by(t1, pa, d).ToPoint(),
            offset_by(t1, pa, -d).ToPoint(),
            offset_by(t2, pa, -d).ToPoint(),
            offset_by(t2, pa, d).ToPoint(),
        });

        auto polygon = std::make_unique<S2Polygon>();
        make_polygon_simple(vertices, *polygon);
        return std::move(polygon);
    }

    double ellipse_rho(const double a, const double b, const double phi)
    {
        return a * b / std::hypot(b * cos(phi), a * sin(phi));
    }

    double ellipse_half_width(const double a,
                              const double b,
                              const double theta,
                              const double position_angle)
    {
        // Use spherical law of sines
        const double rho = ellipse_rho(a, b, theta - position_angle);
        return std::abs(std::asin(sin(rho) * sin(position_angle - theta)));
    }

    vector<S2LatLng> ellipse(const int n,
                             const S2LatLng &center,
                             const double &a,
                             const double &b,
                             const double &theta)
    {
        // all angles in radians
        vector<S2LatLng> el;
        const S1Angle th = S1Angle::Radians(theta);

        assert(n >= 4);
        for (int i = 0; i < n; i++)
        {
            const S1Angle phi = S1Angle::Radians(2 * PI * i / n);
            const S1Angle rho = S1Angle::Radians(ellipse_rho(a, b, phi.radians()));
            el.emplace_back(offset_by(center, th + phi, rho));
        }

        return el;
    }

    std::array<S2LatLng, 2> ellipse_tangent_points(const S2LatLng &center,
                                                   const double a,
                                                   const double b,
                                                   const double theta,
                                                   const double pa)
    {
        // We want the tangent points parallel to the angle pa.
        //
        // Let E be the eccentric anomaly, then the ellipse is:
        //
        //     (x, y) = (a cos(E), b sin(E))
        //
        // Let the angle measured at the center of the ellipse from the
        // semi-major axis to a point on the ellipse be phi, then
        //
        //     tan(E) = a / b * tan(phi)
        //
        // Following
        // https://en.wikipedia.org/wiki/Ellipse#Tangent_slope_as_parameter the
        // point on the ellipse that has a tangent of slope m has:
        //
        //     cot(E) = -m a / b
        //
        // Recognizing that m = tan(pa), then
        //
        //     cot(E) = -(a / b) tan(pa)
        //     --> tan(E) = -(b / a) cot(pa)
        //
        // Now we can relate pa and phi:
        //
        //     tan(phi) = -(b**2 / a**2) cot(pa)
        //
        // But the ellipse is rotated theta degrees, so use pa - theta in
        // place of pa.

        const double phi = std::atan(-b * b / (a * a * tan(pa - theta)));
        const S1Angle rho = S1Angle::Radians(ellipse_rho(a, b, phi));
        const S1Angle xi = S1Angle::Radians(phi + theta);

        // create the tangent points
        const S2LatLng ll1 = offset_by(center, xi, rho);
        const S2LatLng ll2 = offset_by(center, xi, -rho);

        return {ll1, ll2};
    }

    double hav(const double theta)
    {
        return std::pow(sin(theta / 2.0), 2);
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
            A = S1Angle::Radians(PI_2 + cos_c * (PI_2 - position_angle.radians()));

        S1Angle lon = origin.lng() + A;
        S1Angle lat = S1Angle::Radians(std::asin(cos_b));

        return S2LatLng(lat, lon).Normalized();
    }

    double position_angle(const S2Point &a, const S2Point &b)
    {
        return position_angle(S2LatLng(a), S2LatLng(b));
    }

    double position_angle(const S2LatLng &a, const S2LatLng &b)
    {
        if (a.GetDistance(b).radians() < 1e-7)
            return 0;

        // Meeus 1998, Astronomical Algorithms, p116
        S1Angle dra = b.lng() - a.lng();
        return atan2(sin(dra), cos(a.lat()) * tan(b.lat()) - sin(a.lat()) * cos(dra));
    }
}