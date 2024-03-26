#include "config.h"

#include <cstdlib>
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

struct Image
{
    string fov = "";
    bool last_image_of_exposure = false;
    bool last_image_of_night = false;
};

Image new_image()
{
    // repeats a set of observations, drifting at the sidereal rate
    static int tile = 0;      // three tiles per night
    static int iteration = 0; // tile iterations total (used to keep track of night number)
    static int exposure = -1;
    static int image = 0;
    const S1Angle fov_step = S1Angle::Degrees(FOV_WIDTH);
    const S2LatLng image_size = S2LatLng::FromDegrees(IMAGE_WIDTH, IMAGE_WIDTH);
    const S1Angle image_step = S1Angle::Degrees(IMAGE_WIDTH);

    const S1Angle dec0 = S1Angle::Degrees(-30);
    static S1Angle ra0 = S1Angle::Degrees(0);
    static S1Angle dec = dec0;
    static S1Angle ra = ra0;

    const S1Angle two_pi = S1Angle::Radians(2 * PI);
    const S1Angle pi = S1Angle::Radians(PI);

    // note these are the center of the image, not a corner, so there is a 0.5 step offset
    S1Angle ra_i = ra + (image % IMAGES_PER_FOV_WIDTH - IMAGES_PER_FOV_WIDTH / 2.0 + 0.5) * image_step;
    S1Angle dec_i = dec + (image / IMAGES_PER_FOV_WIDTH - IMAGES_PER_FOV_WIDTH / 2.0 + 0.5) * image_step;
    S2LatLngRect fov = S2LatLngRect::FromCenterSize(S2LatLng(dec_i, ra_i).Normalized(), image_size);

    bool last_image_of_exposure = false, last_image_of_night = false;

    image++;
    if (image == IMAGES_PER_FOV_WIDTH * IMAGES_PER_FOV_WIDTH)
    {
        // next image starts a new exposure
        last_image_of_exposure = true;
        exposure++;
        image = 0;

        ra += fov_step / std::cos(dec.radians()); // next RA
        if ((ra - ra0).radians() > RA_COVERAGE)
        {
            // then we're done with this RA stripe
            ra = ra0 + S1Angle::Degrees(235.9 / 86400 * 360 * (iteration / 3));
            dec += fov_step;
        }

        if (dec.degrees() > 90)
        {
            // then we're done with this tile
            iteration++;
            tile = iteration % 3;
            dec = S1Angle::Degrees(-30);
            ra0 = S1Angle::Degrees(235.9 / 86400 * 360 * (iteration / 3));
            ra = ra0;

            if (tile == 0)
                // then we did three tiles, set up a new night
                last_image_of_night = true;
        }
    }

    return Image{format_vertices(fov),
                 last_image_of_exposure,
                 last_image_of_night};
}

void build_test_db()
{
    SBSearch sbs(SBSearch::sqlite3, "sbsearch_test.db", {.log_file = "sbsearch_test.log", .create = true});
    Logger::get_logger().log_level(sbsearch::DEBUG);

    Logger::info() << "Survey setup:"
                   << "\n  Nights: " << NIGHTS
                   << "\n  RA coverage at the equator: " << RA_COVERAGE / DEG << " deg"
                   << "\n  Exposure time: " << EXPOSURE_TIME * 86400 << " s"
                   << "\n  Cadence: " << CADENCE * 86400 << " s"
                   << "\n  Images per exposure: " << IMAGES_PER_FOV_WIDTH * IMAGES_PER_FOV_WIDTH
                   << std::endl;

    Indexer::Options options;
    options.max_spatial_index_cells(MAX_SPATIAL_INDEX_CELLS);
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

    // const double mjd0 = (date_range.first == nullptr) ? 59103.0 : std::ceil(*date_range.second);
    const double mjd0 = 59103.0;
    if (date_range.first == nullptr)
        Logger::info() << "No previous data: starting new survey on mjd = " << mjd0 << std::endl;
    else
    {
        Logger::info() << "Detected prior data: exiting" << std::endl;
        return;
    }

    sbs.drop_observations_indices();

    double mjd;
    int product_id = 0;
    Observations observations;
    observations.reserve(10000);
    sbsearch::ProgressPercent night(NIGHTS);

    while (true) // survey loop
    {
        observations.clear();
        mjd = mjd0 + night.count();
        cout << "night " << night.count() + 1 << std::flush;
        int exposure = 0;

        while (true) // night loop
        {
            Image image = new_image();

            product_id++;

            observations.push_back(
                Observation("test source",
                            "X05",
                            std::to_string(product_id),
                            mjd,
                            mjd + EXPOSURE_TIME,
                            image.fov));

            if (image.last_image_of_night)
                break;

            if (image.last_image_of_exposure)
            {
                exposure++;
                mjd += CADENCE;
            }
        }
        sbs.add_observations(observations);

        night.update();
        night.status();
        if (night.count() == NIGHTS)
            break;
    }

    Logger::info() << "Added " << product_id << " images." << std::endl;
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