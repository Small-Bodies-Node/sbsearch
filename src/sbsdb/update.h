#ifndef SBSDB_UPDATE_H_
#define SBSDB_UPDATE_H_

#include <cinttypes>
#include <vector>

#include "../ephemeris.h"
#include "../found.h"
#include "../observation.h"
#include "../observatory.h"
#include "../moving_target.h"

using std::vector;

namespace sbsearch::sbsdb::update
{
    /**
     * @brief Update a moving target in the database.
     *
     * @param db An sbsearch database instance.
     *
     * @param target The moving target to update.  It must be already be in the
     *               database.
     *
     * Raises MovingTargetError if requirements are not met.
     */
    template <typename DB>
    void moving_target(DB &db, MovingTarget &target);

    /**
     * @brief Update observations in the database.
     *
     * @param db An sbsearch database instance.
     *
     * @param target The observations to update.  `observation_id` must be
     *               defined and in the database.
     *
     * Raises ObservationError if requirements are not met.
     */
    template <typename DB>
    void observations(DB &db, const Observations &obs);

}

#endif // SBSDB_UPDATE_H_
