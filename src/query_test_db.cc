#include "config.h"
#include "test_db.h"

#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <random>

#include "s2/s1angle.h"
#include "s2/s2point.h"
#include "s2/s2latlng.h"

#include "indexer.h"
#include "ephemeris.h"
#include "logging.h"
#include "observation.h"
#include "sbsearch.h"
#include "test_db.h"
#include "util.h"

#define N_COMETS 10

using sbsearch::Ephemeris;
using sbsearch::Found;
using sbsearch::Indexer;
using sbsearch::Observation;
using sbsearch::SBSearch;
using std::cerr;
using std::cout;
using std::endl;
using std::vector;

// ra, dec in deg, rate in deg/day
Ephemeris get_ephemeris(const double mjd0, const double mjd1, const double step,
                        const double ra0, const double dec0,
                        const double ra_rate, const double dec_rate)
{
    static int object_id = 0;
    double ra = ra0, dec = dec0;
    double dec_dir = 1;
    vector<S2Point> vertices;
    vector<double> mjd;

    object_id++;

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

    return Ephemeris(object_id, vertices, mjd, vector<double>(vertices.size(), 1),
                     vector<double>(vertices.size(), 1), vector<double>(vertices.size(), 1));
}

Ephemeris get_random_ephemeris(std::pair<double *, double *> date_range)
{
    // -1000 to 1000 arcsec/hr --> deg/day
    double dec_rate = std::copysign((std::pow(10, 3 * (float)std::rand() / RAND_MAX)) / 3600 * 24, std::rand() - 0.5);
    double ra_rate = std::copysign((std::pow(10, 3 * (float)std::rand() / RAND_MAX)) / 3600 * 24, std::rand() - 0.5);
    double rate = std::hypot(ra_rate, dec_rate);
    printf("- μ = %f deg/day (%f, %f) \n", rate, ra_rate, dec_rate);

    double step = 1.0; // days
    double ra0 = -ra_rate;
    double dec0 = -dec_rate;

    return get_ephemeris(*date_range.first, *date_range.second, step, ra0, dec0, ra_rate, dec_rate);
}

Ephemeris get_fixed_ephemeris(std::pair<double, double> date_range)
{
    double ra_rate, dec_rate;
    ra_rate = FOV_WIDTH / CADENCE;
    dec_rate = 0;
    double rate = std::hypot(ra_rate, dec_rate);
    printf("\n- μ = %f deg/day (%f, %f) \n", rate, ra_rate, dec_rate);

    double step = 1.0; // days
    double ra0 = 0.1 - ra_rate;
    double dec0 = -dec_rate;

    Ephemeris eph = get_ephemeris(date_range.first, date_range.first + 1, step, ra0, dec0, ra_rate, dec_rate);
    return eph;
}

vector<Found> query_sbs(SBSearch *sbs, const Ephemeris &eph)
{
    cout << "  Querying " << eph.num_segments() << " ephemeris segments." << endl;

    auto t0 = std::chrono::steady_clock::now();
    vector<Found> founds = sbs->find_observations(eph);
    std::chrono::duration<double> diff = std::chrono::steady_clock::now() - t0;

    cout << "  Found " << founds.size() << " observation" << (founds.size() == 1 ? "" : "s")
         << " in " << diff.count() << " seconds.\n\n";

    return founds;
}

void query_test_db()
{
    SBSearch sbs(SBSearch::sqlite3, "sbsearch_test.db", "sbsearch_test.log");

    // get date range for query
    auto date_range = sbs.date_range("test source");
    if (date_range.first == nullptr)
        throw std::runtime_error("No observations in database to search.\n");

    cout << "Single point test.\n";
    vector<Observation> observations = sbs.find_observations(S2LatLng::FromDegrees(0, 0.1).ToPoint());
    cout << "\n"
         << observations << "\n";

    cout << "\nGenerating " << (date_range.second - date_range.first) / 365.25 << "-year long ephemerides:\n";

    std::srand(23);
    vector<Found> founds;
    for (int i = 0; i < N_COMETS; ++i)
    {
        Ephemeris eph = get_random_ephemeris(date_range);
        vector<Found> newly_founds = query_sbs(&sbs, eph);
        founds.insert(founds.end(), newly_founds.begin(), newly_founds.end());
    }
    cout << founds;
    cout << "\n\n";
}

int main(int argc, char **argv)
{
    query_test_db();
    return 0;
}