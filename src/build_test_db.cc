#include "config.h"

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>

#include "s2/s1angle.h"
#include "s2/s2latlng.h"
#include "s2/s2latlng_rect.h"

#include "indexer.h"
#include "logging.h"
#include "observation.h"
#include "observatory.h"
#include "sbsearch.h"
#include "test_db.h"
#include "util.h"

using sbsearch::format_vertices;
using sbsearch::Indexer;
using sbsearch::Logger;
using sbsearch::Observation;
using sbsearch::Observations;
using sbsearch::Observatories;
using sbsearch::SBSearch;
using std::cerr;
using std::cout;
using std::endl;
using std::vector;

string new_fov()
{
    // repeats a set of observations, drifting at the sidereal rate
    const int tile_width = int(std::sqrt(TRIPLETS_PER_NIGHT)) + 1;
    static int iteration = 0;
    static int exposure = 0;
    const S2LatLng size = S2LatLng::FromDegrees(FOV_WIDTH, FOV_WIDTH);
    const S1Angle step = S1Angle::Degrees(FOV_WIDTH);
    const S1Angle two_pi = S1Angle::Radians(2 * PI);
    const S1Angle pi = S1Angle::Radians(PI);

    exposure++;
    int tile_exposure = exposure % TRIPLETS_PER_NIGHT;
    if (tile_exposure == 0)
        iteration++;

    S1Angle ra = (tile_exposure % tile_width - 1) * step + S1Angle::Degrees(235.9 / 86400 * 360 * (iteration / 3));
    S1Angle dec = (tile_exposure / tile_width) * step;

    S2LatLngRect rect = S2LatLngRect::FromCenterSize(S2LatLng(dec, ra).Normalized(), size);
    return format_vertices(rect);
}

void build_test_db()
{
    SBSearch sbs(SBSearch::sqlite3, "sbsearch_test.db", {.log_file = "sbsearch_test.log"});

    Logger::info() << "Survey setup:"
                   << "\n  Nights: " << NIGHTS
                   << "\n  Exposures per night: " << EXPOSURES_PER_NIGHT
                   << "\n  Cadence: " << CADENCE * 86400 << " s"
                   << "\n  Exposure time: " << EXPOSURE * 86400 << " s" << std::endl;

    Indexer::Options options;
    options.max_spatial_cells(MAX_SPATIAL_CELLS);
    options.max_spatial_resolution(MAX_SPATIAL_RESOLUTION);
    options.min_spatial_resolution(MIN_SPATIAL_RESOLUTION);
    options.temporal_resolution(TEMPORAL_RESOLUTION);

    auto date_range = sbs.db()->observation_date_range();

    Indexer::Options options_saved = sbs.indexer_options();

    // make sure database options match what we think they should be
    if (options != options_saved)
    {
        // If they do not match and there are observations in the database, throw an error.
        if (date_range.first != nullptr)
            throw std::runtime_error("Configuration does not match database: re-index before adding more data.");

        // otherwise, quietly update them
        Logger::debug() << "Updating database configuration." << std::endl;
        sbs.reindex(options);
    }

    // and add our observatory
    Observatories observatories = sbs.db()->get_observatories();
    if (observatories.find("X05") == observatories.end())
        sbs.db()->add_observatory("X05", {289.25058, 0.864981, -0.500958});

    const double mjd0 = (date_range.first == nullptr) ? 59103.0 : std::ceil(*date_range.second);
    if (date_range.first == nullptr)
        Logger::info() << "No previous data: starting new survey on mjd = " << mjd0 << std::endl;
    else
        Logger::info() << "Detected prior data: appending observations and starting with mjd = " << mjd0 << std::endl;

    sbs.drop_observations_indices();

    double mjd;
    int product_id = 0;
    Observations observations;
    observations.reserve(EXPOSURES_PER_NIGHT);
    sbsearch::ProgressPercent progress(EXPOSURES_PER_NIGHT * NIGHTS);
    for (int night = 0; night < NIGHTS; night++)
    {
        observations.clear();
        mjd = mjd0 + night;
        cout << "\rnight " << night + 1 << std::flush;
        for (int triplet = 0; triplet < 3; triplet++)
        {
            for (int i = 0; i < TRIPLETS_PER_NIGHT; ++i)
            {
                product_id++;
                observations.push_back(Observation("test source",
                                                   "X05",
                                                   std::to_string(product_id),
                                                   mjd + CADENCE * i,
                                                   mjd + CADENCE * i + EXPOSURE,
                                                   new_fov()));
            }
            mjd += CADENCE * (1 + TRIPLETS_PER_NIGHT * (triplet + 1));
        }
        sbs.add_observations(observations);
        progress.update(observations.size());
        progress.status();
    }
    Logger::info() << "Added " << progress.count() << " observations" << std::endl;

    Logger::debug() << "Creating indices." << std::endl;
    sbs.create_observations_indices();
    Logger::info() << "Done." << std::endl;
}

int main(int argc, char **argv)
{
    try
    {
        build_test_db();
    }
    catch (const std::runtime_error &error)
    {
        Logger::error() << error.what() << std::endl;
    }
    return 0;
}