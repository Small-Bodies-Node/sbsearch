#include "config.h"

#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <random>
#include "s2/s1angle.h"
#include "s2/s2point.h"
#include "s2/s2latlng.h"

#include "date.h"
#include "horizons.h"
#include "indexer.h"
#include "ephemeris.h"
#include "found.h"
#include "logging.h"
#include "moving_target.h"
#include "observation.h"
#include "sbsearch.h"
#include "test_db.h"
#include "util.h"

#define N_COMETS 30
#define FIXED_SEARCH false
#define ENCKE_SEARCH false
#define PARALLAX_SEARCH true
#define MAX_SPATIAL_QUERY_CELLS 8

using namespace sbsearch;
using std::cerr;
using std::cout;
using std::endl;
using std::vector;

// ra, dec in deg, rate in deg/day
Ephemeris get_ephemeris(const double mjd0, const double mjd1, const double step,
                        const double ra0, const double dec0,
                        const double ra_rate, const double dec_rate,
                        const double delta, const double tp)
{
    static int moving_target_id = 0;
    double ra = ra0, dec = dec0;
    double dec_dir = 1;
    Ephemeris::Data data;

    moving_target_id++;

    for (double mjd = mjd0 - step; mjd < mjd1 + step; mjd += step)
    {
        data.push_back({
            .mjd = mjd,
            .tmtp = mjd - tp,
            .ra = ra,
            .dec = dec,
            .rh = delta + 1,
            .delta = delta,
            .phase = 1,
        });
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

    MovingTarget target("Target " + std::to_string(moving_target_id), moving_target_id);
    return Ephemeris(target, data);
}

Ephemeris get_random_ephemeris(std::pair<double *, double *> date_range)
{
    // -1000 to 1000 arcsec/hr --> deg/day
    const double dec_rate = std::copysign((std::pow(10, 3 * (float)std::rand() / RAND_MAX)) / 3600 * 24, std::rand() - 0.5);
    const double ra_rate = std::copysign((std::pow(10, 3 * (float)std::rand() / RAND_MAX)) / 3600 * 24, std::rand() - 0.5);
    const double rate = std::hypot(ra_rate, dec_rate);
    const double delta = std::pow(10, (float)std::rand() / RAND_MAX - 1);
    const double tp = ((float)std::rand() / RAND_MAX - 0.5) * 365 * 5 + 59103.0;
    // printf("- μ = %f deg/day (%f, %f)\n", rate, ra_rate, dec_rate);
    // printf("- Delta = %f au\n", delta);
    Logger::info() << "- μ = " << rate << " deg/day (" << ra_rate << ", " << dec_rate << ")" << endl;
    Logger::info() << "- Δ = " << delta << " au" << endl;

    const double step = 1; // days
    const double ra0 = -ra_rate;
    const double dec0 = -dec_rate;

    return get_ephemeris(*date_range.first, *date_range.second, step, ra0, dec0, ra_rate, dec_rate, delta, tp);
}

Ephemeris get_fixed_ephemeris(std::pair<double *, double *> date_range)
{
    double ra_rate, dec_rate;
    ra_rate = FOV_WIDTH / CADENCE;
    dec_rate = 0;
    double rate = std::hypot(ra_rate, dec_rate);
    printf("\n- μ = %f deg/day (%f, %f) \n", rate, ra_rate, dec_rate);

    double step = 1.0; // days
    double ra0 = 0.1 - ra_rate;
    double dec0 = -dec_rate;

    Ephemeris eph = get_ephemeris(*date_range.first, *date_range.second, step, ra0, dec0, ra_rate, dec_rate, 1, 59103.0);
    return eph;
}

vector<Found> query_sbs(SBSearch *sbs, const Ephemeris &eph)
{
    Logger::info() << "Querying " << eph.num_segments() << " ephemeris segments." << endl;

    auto t0 = std::chrono::steady_clock::now();
    vector<Found> founds = sbs->find_observations(eph, {.parallax = PARALLAX_SEARCH, .max_spatial_query_cells = MAX_SPATIAL_QUERY_CELLS});
    std::chrono::duration<double> diff = std::chrono::steady_clock::now() - t0;

    Logger::info() << "Found " << founds.size() << " observation" << (founds.size() == 1 ? "" : "s")
                   << " in " << diff.count() << " seconds." << endl;

    return founds;
}

void query_test_db()
{
    SBSearch sbs(SBSearch::sqlite3, "sbsearch_test.db", {.log_file = "sbsearch_test.log", .log_level = LogLevel::DEBUG});

    // get date range for query
    const std::pair<double *, double *>
        date_range = sbs.db()->observation_date_range("test source");
    if (date_range.first == nullptr)
        throw std::runtime_error("No observations in database to search.\n");

    if (FIXED_SEARCH)
    {
        Logger::info() << "Single point test." << endl;
        Observations observations = sbs.find_observations(S2LatLng::FromDegrees(0, 0.1).ToPoint());
        observations[0].format.show_fov = true;
        cout << "\n"
             << observations << "\n";
    }

    Logger::info() << "Generating " << (*date_range.second - *date_range.first) / 365.25
                   << "-year long ephemerides" << endl;

    std::srand(23);
    vector<Found> founds;
    for (int i = 0; i < N_COMETS; ++i)
    {
        Ephemeris eph = get_random_ephemeris(date_range);
        vector<Found> newly_founds = query_sbs(&sbs, eph);
        founds.insert(founds.end(), newly_founds.begin(), newly_founds.end());
    }

    if (ENCKE_SEARCH)
    {
        MovingTarget encke("2P");
        Horizons horizons(encke,
                          "500@399",
                          Date(*date_range.first),
                          Date(*date_range.second),
                          "1d",
                          true);
        Ephemeris eph(encke, horizons.get_ephemeris_data());
        vector<Found> newly_founds = query_sbs(&sbs, eph);
        founds.insert(founds.end(), newly_founds.begin(), newly_founds.end());
    }
    cout << "\n";
    cout << founds;
    cout << "\n\n";
}

int main(int argc, char **argv)
{
    query_test_db();
    return 0;
}