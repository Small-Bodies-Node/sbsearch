#ifndef TEST_DB_H_
#define TEST_DB_H_

#include "util.h"

#define FOV_WIDTH 7.0 /* degrees */
#define IMAGES_PER_FOV_WIDTH 4
#define IMAGE_WIDTH (FOV_WIDTH / IMAGES_PER_FOV_WIDTH)

#define NIGHTS 365
// RA coverage on the equator, may not work if RA_COVERAGE + FOV_WIDTH > 180 deg
#define RA_COVERAGE (160 * DEG)
#define IMAGES_PER_EXPOSURE (IMAGES_PER_FOV_WIDTH * IMAGES_PER_FOV_WIDTH)

#define EXPOSURE_TIME (30.0 / 86400)
#define CADENCE (45.0 / 86400)

#define MAX_SPATIAL_CELLS 8
#define MAX_SPATIAL_RESOLUTION (1 * DEG)
#define MIN_SPATIAL_RESOLUTION (5 * ARCMIN)
#define TEMPORAL_RESOLUTION 1 /* per day */

#endif // TEST_DB_H_