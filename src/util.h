#ifndef SBSEARCH_UTIL_H_
#define SBSEARCH_UTIL_H_

#include <vector>
#include <set>
#include <string>
#include <s2/s2latlng.h>
#include <s2/s2latlng_rect.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>

#define PI 3.14159265358979323846
#define PI_2 1.57079632679489661923
#define DEG (PI / 180)
#define ARCMIN (PI / 10800)
#define ARCSEC (PI / 648000)
#define HOUR (PI / 12)

#define UNDEF_TIME -1
#define UNDEF_ANGLE -999
#define UNDEF_UNC -1

#define CERR(x) (std::cerr << x << std::endl)

using std::string;
using std::vector;

namespace sbsearch
{
    // calculate the position angle (angle from North) from point a to point b
    double position_angle(const S2Point &a, const S2Point &b);

    // offset `distance` from `point` along `position_angle`.
    S2LatLng offset_by(const S2LatLng &point, const S1Angle &position_angle, const S1Angle &distance);

    // Generate an ellipse, composed of n points, angles in radians.
    vector<S2LatLng> ellipse(const int n, const S2LatLng &center, const double &a, const double &b, const double &theta);

    // Split a string given delimiter.  The delimiter is not included in the
    // output.
    vector<string> split(string s, const char delimiter);

    // Join a vector of strings with the delimiter.
    string join(const vector<string> s, const char *delimiter);

    // vector values must always be increasing
    bool is_increasing(const vector<double> &v);

    // String-formatted vertices, comma-separated RA:Dec pairs in units of
    // degrees, e.g., "0:0, 0:1, 1:1".  For a polygon, only the first loop is
    // checked.
    string format_vertices(const vector<S2LatLng> vertices);
    string format_vertices(const vector<S2Point> vertices);
    string format_vertices(const S2LatLngRect fov);
    string format_vertices(const S2Polygon &polygon);

    // RA, dec in units of degrees
    string format_vertices(int num_vertices, const double *ra, const double *dec);

    // Convert string format ("RA:Dec, ...", units of degrees) to vector of points
    vector<S2Point> make_vertices(string str);

    void make_polygon(const vector<S2Point> &vertices, S2Polygon &polygon);
    void make_polygon(string str, S2Polygon &polygon);

    // Add a padding around the polygon given by pad in arcmin.  Padding must be
    // > 0 or else the original polygon will be returned unmodified.
    void padded_polygon(const S2Polygon &polygon, const double padding, S2Polygon &result);

    // interpolation
    // there are no limits on `frac`, use `frac` < 0 or > 1 to extrapolate.
    inline double interp(const double a, const double b, const double frac)
    {
        return a + (b - a) * frac;
    }

    string mjd2cal(const double &mjd);
}

#endif // SBSEARCH_UTIL_H_