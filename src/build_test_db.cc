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
#include "sbsdb/add.h"
#include "sbsdb/get.h"
#include "sbsdb/postgresql.h"
#include "test_db.h"
#include "util/string.h"

using namespace sbsearch;
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

    return Image{sbsearch::util::format_vertices(fov),
                 last_image_of_exposure,
                 last_image_of_night};
}

template <typename DB>
void build_test_db(string url)
{
    SBSearch<DB> sbs(url, {.log_file = "sbsearch_test.log", .create = true});
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

    const auto date_range = sbsdb::get::observations_date_range(sbs.db());

    Indexer::Options options_saved = sbs.indexer_options();

    // make sure database options match what we think they should be
    if (options != options_saved)
    {
        cout << "\nCurrent index setup:"
             << "\n  Maximum spatial cells: " << options_saved.max_spatial_index_cells()
             << "\n  Minimum spatial level: "
             << options_saved.min_spatial_level()
             << " (" << options_saved.max_spatial_resolution() / DEG << " deg)"
             << "\n  Maximum spatial level: "
             << options_saved.max_spatial_level()
             << " (" << options_saved.min_spatial_resolution() / DEG << " deg)"
             << "\n  Temporal resolution (1/day): " << options_saved.temporal_resolution()
             << "\n\nExpected index setup:"
             << "\n  Maximum spatial index cells: " << options.max_spatial_index_cells()
             << "\n  Minimum spatial level: "
             << options.min_spatial_level()
             << " (" << options.max_spatial_resolution() / DEG << " deg)"
             << "\n  Maximum spatial level: "
             << options.max_spatial_level()
             << " (" << options.min_spatial_resolution() / DEG << " deg)"
             << "\n  Temporal resolution (1/day): " << options.temporal_resolution()
             << "\n\n";

        // If they do not match and there are observations in the database, throw an error.
        if (date_range.first)
            throw std::runtime_error("Configuration does not match database: re-index before adding more data.");

        // otherwise, quietly update them
        Logger::debug() << "Updating database configuration." << std::endl;
        sbs.reindex(options);
    }

    // and add our observatory
    Observatories observatories = sbsdb::get::all_observatories(sbs.db());
    if (observatories.find("X05") == observatories.end())
        sbsdb::add::observatory(sbs.db(), {289.25058, 0.864981, -0.500958, "X05"});

    // const double mjd0 = (!date_range.first) ? 59103.0 : std::ceil(date_range.second);
    const double mjd0 = 59103.0;
    if (!date_range.first)
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
    observations.data.reserve(10000);
    sbsearch::ProgressPercent night(NIGHTS);

    while (true) // survey loop
    {
        observations.data.clear();
        mjd = mjd0 + night.count();
        int exposure = 0;

        while (true) // night loop
        {
            Image image = new_image();

            product_id++;

            observations.append(
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
        build_test_db<sbsdb::Postgresql>("postgres:///sbsearch_test_db");
    }
    catch (const std::runtime_error &error)
    {
        Logger::error() << error.what() << std::endl;
    }
    return 0;
}