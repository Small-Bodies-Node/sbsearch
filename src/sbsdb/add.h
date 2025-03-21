#ifndef SBSDB_ADD_H_
#define SBSDB_ADD_H_

#include <cinttypes>
#include <vector>

#include "../ephemeris.h"
#include "../found.h"
#include "../observation.h"
#include "../observatory.h"
#include "../moving_target.h"

using std::vector;

namespace sbsearch::sbsdb::add
{
    /**
     * @brief Add a moving target to the database.
     *
     * @param db An sbsearch database instance.
     *
     * @param target The moving target to add, `moving_target_id` must be null
     *               and will be updated with the new ID.
     */
    template <typename DB>
    void moving_target(DB &db, MovingTarget &target);

    /**
     * @brief Add observations to the database.
     *
     * @param db An sbsearch database instance.
     *
     * @param target The observations to add.  `observation_id` must be null and
     *               `terms` must be defined.
     *
     * @returns Observations A copy of the observations with updated
     *                       `observation_id`.
     *
     * Raises ObservationError if observation requirements are not met.
     */
    template <typename DB>
    Observations observations(DB &db, const Observations &obs);

}

#endif // SBSDB_ADD_H_
