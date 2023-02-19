#ifndef TEST_DB_H_
#define TEST_DB_H_

#include "util.h"

#define FOV_WIDTH 0.5 /* degrees */
#define NIGHTS 365
#define TRIPLETS_PER_NIGHT 333
#define EXPOSURES_PER_NIGHT (TRIPLETS_PER_NIGHT * 3)
#define CADENCE (0.5 / EXPOSURES_PER_NIGHT)
#define EXPOSURE (0.8 * CADENCE)

#define MAX_SPATIAL_CELLS 8
#define MAX_SPATIAL_RESOLUTION (1 * DEG)
#define MIN_SPATIAL_RESOLUTION (1 * ARCMIN)
#define TEMPORAL_RESOLUTION 10 /* per day */

#endif // TEST_DB_H_