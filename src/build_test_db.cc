#include "config.h"

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>

#include "s2/s1angle.h"
#include "s2/s2latlng.h"
#include "s2/s2latlng_rect.h"

#include "indexer.h"
#include "observation.h"
#include "sbsearch.h"
#include "test_db.h"
#include "util.h"

using sbsearch::format_vertices;
using sbsearch::Indexer;
using sbsearch::Observation;
using sbsearch::SBSearch;
using std::cerr;
using std::cout;
using std::endl;
using std::vector;

string new_fov()
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

    S2LatLngRect rect = S2LatLngRect::FromCenterSize(S2LatLng(dec, ra).Normalized(), size);
    return format_vertices(rect);
}

void build_test_db()
{
    Indexer::Options options;
    options.max_spatial_cells(MAX_SPATIAL_CELLS);
    options.max_spatial_resolution(MAX_SPATIAL_RESOLUTION);
    options.min_spatial_resolution(MIN_SPATIAL_RESOLUTION);
    options.temporal_resolution(TEMPORAL_RESOLUTION);
    SBSearch sbs(SBSearch::sqlite3, "sbsearch_test.db", options);

    cout << "\nSurvey setup:"
         << "\n  Nights: " << NIGHTS
         << "\n  Exposures per night: " << EXPOSURES_PER_NIGHT
         << "\n  Cadence: " << CADENCE * 86400 << " s"
         << "\n  Exposure time: " << EXPOSURE * 86400 << " s\n\n";

    sbs.drop_observations_indices();

    const double mjd0 = 59103.0;
    double mjd;
    int product_id = 0;
    vector<Observation> observations;
    observations.reserve(EXPOSURES_PER_NIGHT);
    for (int night = 0; night < NIGHTS; night++)
    {
        observations.clear();
        mjd = mjd0 + night;
        cout << "\rnight " << night + 1 << std::flush;
        for (int i = 0; i < EXPOSURES_PER_NIGHT; ++i)
        {
            mjd += CADENCE;
            product_id++;
            observations.push_back(Observation("test source", std::to_string(product_id), mjd, mjd + EXPOSURE, new_fov()));
        }
        sbs.add_observations(observations);
    }

    cout << endl
         << endl;

    cout << "Creating indices.\n";
    sbs.create_observations_indices();
}

int main(int argc, char **argv)
{
    build_test_db();
    return 0;
}