#ifndef DB3_H_
#define DB3_H_

#include <vector>
#include <string>
#include <sqlite3.h>

#include <s2/s2point.h>
#include <s2/s2polygon.h>

#define PI 3.14159265358979323846
#define PI_2 1.57079632679489661923
#define DEG PI / 180
#define ARCMIN PI / 10800
#define ARCSEC PI / 648000

#define CERR(x) (std::cerr << x << std::endl)

using std::string;
using std::vector;

namespace sbsearch
{
    void sql_check(int rc, char *error_message);
    void sql_execute(sqlite3 *db, const char *statement);
    vector<string> mjd_to_time_terms(const double start, const double stop);
    double position_angle(const S2Point &a, const S2Point &b);
    vector<string> split(string s, const char delimiter);
    string join(const vector<string> s, const char *delimiter);
    vector<S2Point> makeVertices(string str);
    std::unique_ptr<S2Polygon> makePolygon(vector<S2Point> vertices);
    std::unique_ptr<S2Polygon> makePolygon(string str);
}
#endif