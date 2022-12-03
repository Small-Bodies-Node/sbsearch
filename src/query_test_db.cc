#include "test_db.h"
#include "sbsearch.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <sstream>
#include <iostream>
#include <list>
#include <vector>
#include <string>
#include <random>
#include <ctime>

#include <sqlite3.h>

#include "s2/s2region_term_indexer.h"
#include "s2/s2polyline.h"
#include "s2/s2metrics.h"
#include "s2/s2latlng.h"
#include "s2/s2builder.h"
#include "s2/s2error.h"
#include "s2/s2testing.h"
#include "s2/s2builderutil_s2polygon_layer.h"
#include "s2/s2text_format.h"

#include "util.h"
#include "observation.h"
#include "ephemeris.h"
#include "sbsdb_sqlite3.h"

#define N_COMETS 10

using sbsearch::Ephemeris;
using sbsearch::Found;
using sbsearch::Observation;
using sbsearch::SBSearchDatabaseSqlite3;
using std::cerr;
using std::cout;
using std::endl;
using std::vector;

Ephemeris get_ephemeris(const double mjd0, const double mjd1, const double step,
                        const double ra0, const double dec0,
                        const double ra_rate, const double dec_rate)
{
    double ra = ra0, dec = dec0;
    double dec_dir = 1;
    vector<S2Point> all_vertices;
    vector<double> all_times;
    const S1Angle tol = S1Angle::Radians(EPHEMERIS_TOLERANCE);

    for (double mjd = mjd0; mjd <= mjd1; mjd += step)
    {
        all_vertices.push_back(S2LatLng::FromDegrees(dec, ra).ToPoint());
        all_times.push_back(mjd);
        ra += ra_rate * std::cos(dec * PI / 180) * step;
        dec += dec_dir * dec_rate * step;
        if (dec > 90)
        {
            dec = 90 - (dec - 90);
            dec_dir *= -1;
            ra += 180;
        }
        else if (dec < -90)
        {
            dec = -90 - (dec + 90);
            dec_dir *= -1;
            ra += 180;
        }

        if (ra > 180)
        {
            ra -= 360;
        }
        else if (ra < -180)
        {
            ra += 360;
        }
    }
    S2Polyline full_resolution(all_vertices);

    // simplify
    vector<int> indices;
    vector<S2Point> vertices;
    vector<double> times;
    full_resolution.SubsampleVertices(tol, &indices);
    for (const auto &i : indices)
    {
        vertices.push_back(all_vertices[i]);
        times.push_back(all_times[i]);
    }

    return Ephemeris(vertices, times);
}

Ephemeris get_random_ephemeris(std::pair<double, double> date_range)
{
    double step = 1.0; // days
    double ra0 = FOV_WIDTH * 1.5;
    double dec0 = 0;
    double ra_rate, dec_rate;
    // -100 to 100 arcsec/hr
    dec_rate = ((float)std::rand() / RAND_MAX - 0.5) * 200 / 3600 * 24; // deg/day
    ra_rate = ((float)std::rand() / RAND_MAX - 0.5) * 200 / 3600 * 24;
    printf("- μ = %lf deg/day", std::hypot(ra_rate, dec_rate));

    return get_ephemeris(date_range.first, date_range.second, step, ra0, dec0, ra_rate, dec_rate);
}

std::pair<double, double> survey_date_range(SBSearchDatabaseSqlite3 &db)
{
    double mjd_start, mjd_stop;
    mjd_start = db.get_one_value<double>("SELECT MIN(mjd_start) FROM observations;");
    mjd_stop = db.get_one_value<double>("SELECT MAX(mjd_stop) FROM observations;");

    return std::pair<double, double>(mjd_start, mjd_stop);
}

void query_test_db()
{
    // get date range for query
    SBSearchDatabaseSqlite3 db("sbsearch_test.db");
    std::pair<double, double> date_range = survey_date_range(db);
    cout << "\nGenerating " << (date_range.second - date_range.first) / 365.25 << " year long ephemerides:\n";

    // std::srand((unsigned)std::time(0));
    std::srand(23);
    for (int i = 0; i < N_COMETS; ++i)
    {
        Ephemeris eph = get_random_ephemeris(date_range);
        cout << "\n  Querying " << eph.num_segments() << " ephemeris segments." << endl;
        vector<Found> found = db.find_observations(eph);
        cout << "  Found " << found.size() << " observations." << endl;
        std::for_each(found.begin(), found.begin() + std::min<int>(found.size(), 30), [](auto &f)
                      { cout << f.obs.observation_id() << " " << f.obs.mjd_start() << " " << f.obs.mjd_stop() << " " << f.obs.fov() << endl; });
        if (found.size() > 30)
            cout << "...\n";
    }
    cout << "\n\n";
}

int main(int argc, char **argv)
{
    query_test_db();
    return 0;
}