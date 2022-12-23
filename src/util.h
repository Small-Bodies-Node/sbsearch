#ifndef SBSEARCH_UTIL_H_
#define SBSEARCH_UTIL_H_

#include <vector>
#include <string>
#include <sqlite3.h>

#include <s2/s2latlng.h>
#include <s2/s2latlng_rect.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>

#define PI 3.14159265358979323846
#define PI_2 1.57079632679489661923
#define DEG (PI / 180)
#define ARCMIN (PI / 10800)
#define ARCSEC (PI / 648000)

#define CERR(x) (std::cerr << x << std::endl)

using std::string;
using std::vector;

namespace sbsearch
{
    // calculate the position angle (angle from North) from point a to point b
    double position_angle(const S2Point &a, const S2Point &b);

    // offset `distance` from `point` along `position_angle`.
    S2LatLng offset_by(const S2LatLng &point, const S1Angle &position_angle, const S1Angle &distance);

    // generate an ellipse, composed of n points
    vector<S2LatLng> ellipse(const int n, const S2LatLng &center, const double &a, const double &b, const double &theta);

    // split and join vectors of strings
    vector<string> split(string s, const char delimiter);
    string join(const vector<string> s, const char *delimiter);

    // string-formatted vertices
    // string format is comma-separated RA:Dec pairs in units of degrees, e.g., "0:0, 0:1, 1:1"
    string format_vertices(vector<S2LatLng> vertices);
    string format_vertices(vector<S2Point> vertices);
    string format_vertices(S2LatLngRect fov);
    // units of degrees
    string format_vertices(int num_vertices, double *ra, double *dec);

    // Convert string format ("RA:Dec, ...", units of degrees) to vector of points
    vector<S2Point> makeVertices(string str);

    void makePolygon(const vector<S2Point> &vertices, S2Polygon &polygon);
    void makePolygon(string str, S2Polygon &polygon);

    // interpolation
    // there are no limits on `frac`, use `frac` < 0 or > 1 to extrapolate.
    inline double interp(const double a, const double b, const double frac)
    {
        return a + (b - a) * frac;
    }
}
#endif // SBSEARCH_UTIL_H_