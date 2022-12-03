#include "test_db.h"
#include "sbsearch.h"

#include <cstdio>
#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include <numeric>

#include <sqlite3.h>

#include "s2/s2region_term_indexer.h"
#include "s2/s2latlng.h"
#include "s2/s2latlng_rect.h"
#include "s2/s2metrics.h"

#include "observation.h"
#include "util.h"
#include "sbsdb_sqlite3.h"

using sbsearch::mjd_to_time_terms;
using sbsearch::Observation;
using sbsearch::SBSearchDatabaseSqlite3;
using std::cerr;
using std::cout;
using std::endl;
using std::vector;

void new_fov(S2LatLngRect &fov)
{
    const static S2LatLng size = S2LatLng::FromDegrees(FOV_WIDTH, FOV_WIDTH);
    const static S1Angle step = S1Angle::Degrees(FOV_WIDTH * 2);
    const static S1Angle two_pi = S1Angle::Radians(2 * PI);
    const static S1Angle pi = S1Angle::Radians(PI);
    static S1Angle ra = S1Angle::Zero();
    static S1Angle dec = S1Angle::Zero();

    ra += step;
    if (ra.radians() > PI)
    {
        dec += step;
        ra -= two_pi;
    }

    if (dec.radians() > PI_2)
    {
        dec -= pi;
    }

    fov = S2LatLngRect::FromCenterSize(S2LatLng(dec, ra).Normalized(), size);
}

void build_test_db()
{
    SBSearchDatabaseSqlite3 db("sbsearch_test.db");
    db.setup_tables();
    cout << "\nSurvey setup:"
         << "\n  Nights: " << NIGHTS
         << "\n  Exposures per night: " << EXPOSURES_PER_NIGHT
         << "\n  Cadence: " << CADENCE * 86400 << " s"
         << "\n  Exposure time: " << EXPOSURE * 86400 << " s\n\n";

    db.drop_time_indices();

    const double mjd0 = 59103.0;
    double mjd;
    vector<Observation> observations;
    observations.reserve(EXPOSURES_PER_NIGHT);
    S2LatLngRect fov;
    for (int night = 0; night < NIGHTS; night++)
    {
        observations.clear();
        mjd = mjd0 + night;
        cout << "\rnight " << night + 1 << std::flush;
        for (int i = 0; i < EXPOSURES_PER_NIGHT; ++i)
        {
            mjd += CADENCE;
            new_fov(fov);
            observations.push_back(Observation(mjd, mjd + EXPOSURE, fov));
        }
        db.add_observations(observations);
    }

    cout << endl
         << endl;

    cout << "Creating indices.\n";
    db.setup_tables();
}

int main(int argc, char **argv)
{
    build_test_db();
    return 0;
}