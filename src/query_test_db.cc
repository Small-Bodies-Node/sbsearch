#include "config.h"
#include "test_db.h"

#include <cmath>
#include <ctime>
#include <iostream>
#include <vector>
#include <string>
#include <random>
// #include <algorithm>
// #include <cstdio>
// #include <cstdlib>
// #include <set>
// #include <sstream>
// #include <list>

// #include <sqlite3.h>

#include "s2/s1angle.h"
#include "s2/s2point.h"
// #include "s2/s2region_term_indexer.h"
// #include "s2/s2polyline.h"
// #include "s2/s2metrics.h"
// #include "s2/s2latlng.h"
// #include "s2/s2builder.h"
// #include "s2/s2error.h"
// #include "s2/s2testing.h"
// #include "s2/s2builderutil_s2polygon_layer.h"
// #include "s2/s2text_format.h"

#include "indexer.h"
#include "ephemeris.h"
#include "observation.h"
#include "sbsearch.h"
#include "test_db.h"
#include "util.h"

#define N_COMETS 1

using sbsearch::Ephemeris;
using sbsearch::Found;
using sbsearch::Indexer;
using sbsearch::Observation;
using sbsearch::SBSearch;
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
    vector<S2Point> vertices;
    vector<double> mjd;
    const S1Angle tol = S1Angle::Radians(1 * ARCSEC);

    for (double t = mjd0 - step; t < mjd1 + step; t += step)
    {
        vertices.push_back(S2LatLng::FromDegrees(dec, ra).ToPoint());
        mjd.push_back(t);
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

    return Ephemeris(vertices, mjd, vector<double>(vertices.size(), 1),
                     vector<double>(vertices.size(), 1), vector<double>(vertices.size(), 1));
}

Ephemeris get_random_ephemeris(std::pair<double, double> date_range)
{
    double step = 1.0; // days
    double ra0 = FOV_WIDTH * 1.5;
    double dec0 = 0;
    double ra_rate, dec_rate;
    // -100 to 100 arcsec/hr
    // dec_rate = ((float)std::rand() / RAND_MAX - 0.5) * 200 / 3600 * 24; // deg/day
    // ra_rate = ((float)std::rand() / RAND_MAX - 0.5) * 200 / 3600 * 24;
    ra_rate = FOV_WIDTH / CADENCE;
    dec_rate = 0;
    printf("- μ = %f deg/day (%f, %f) \n", std::hypot(ra_rate, dec_rate), ra_rate, dec_rate);

    return get_ephemeris(date_range.first, date_range.second, step, ra0, dec0, ra_rate, dec_rate);
}

void print_found(Found found)
{
    printf("\"%s\" \"%s\" %.6f \"%s\" %.6f %.6f %.3f %.3f %.2f\n",
           found.observation.source().c_str(), found.observation.product_id().c_str(), found.observation.mjd_start(),
           found.observation.fov().c_str(), found.ephemeris.ra(0), found.ephemeris.dec(0), found.ephemeris.rh(0),
           found.ephemeris.delta(0), found.ephemeris.phase(0));
}

void print_ephemeris(Ephemeris eph)
{
    for (int i = 0; i < eph.num_vertices(); i++)
    {
        printf("%.6f %.6f %.6f %.3f %.3f %.2f\n",
               eph.mjd(i), eph.ra(i), eph.dec(i), eph.rh(i), eph.delta(i),
               eph.phase(i));
    }
}

void query_test_db()
{
    Indexer::Options options;
    options.max_spatial_cells(MAX_SPATIAL_CELLS);
    options.max_spatial_resolution(MAX_SPATIAL_RESOLUTION);
    options.min_spatial_resolution(MIN_SPATIAL_RESOLUTION);
    SBSearch sbs(SBSearch::sqlite3, "sbsearch_test.db", options);

    // get date range for query
    std::pair<double, double> date_range = sbs.date_range("test source");
    cout << "\nGenerating " << (date_range.second - date_range.first) / 365.25 << " year long ephemerides:\n";

    // std::srand(23);
    for (int i = 0; i < N_COMETS; ++i)
    {
        Ephemeris eph = get_random_ephemeris(date_range);
        print_ephemeris(eph);
        cout << "\n  Querying " << eph.num_segments() << " ephemeris segments." << endl;
        vector<Found> found = sbs.find_observations(eph);
        cout << "  Found " << found.size() << " observations." << endl;
        std::for_each(found.begin(), found.begin() + std::min<int>(found.size(), 30), print_found);

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