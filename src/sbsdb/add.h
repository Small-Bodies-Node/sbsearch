#ifndef SBSDB_ADD_H_
#define SBSDB_ADD_H_

#include <cinttypes>
#include <vector>

#include "found.h"
#include "observation.h"
#include "observatory.h"
#include "moving_target.h"
#include "ephemeris/ephemeris.h"

using sbsearch::ephemeris::Ephemeris;
using std::vector;

namespace sbsearch::sbsdb::add
{
    /**
     * @brief Add an ephemeris to the database.
     *
     * The target must be in the database.
     *
     * @param db An sbsearch database instance.
     *
     * @param eph The ephemeris to add.
     *
     * @throw MovingTargetError when the database target does not match the ephemeris target.
     */
    template <typename DB>
    void ephemeris(DB *db, Ephemeris &eph);

    /**
     * @brief Add an entry to the found object table.
     *
     * The target must be in the database.
     *
     * @param db An sbsearch database instance.
     *
     * @param founds The found data to add.
     *
     * @throw MovingTargetError when the database target does not match the ephemeris target.
     */
    template <typename DB>
    void found(DB *db, const Founds &founds);

    /**
     * @brief Add a moving target to the database.
     *
     * @param db An sbsearch database instance.
     *
     * @param target The moving target to add, `moving_target_id` must be null
     *               and will be updated with the new ID.
     */
    template <typename DB>
    void moving_target(DB *db, MovingTarget &target);

    /**
     * @brief Add observations to the database.
     *
     * Switches to many_observations when there are many observations.
     *
     * @param db An sbsearch database instance.
     *
     * @param obs The observations to add.  `observation_id` must be null and
     *            `terms` must be defined.  `observation_id` will be updated.
     *
     * Raises ObservationError if observation requirements are not met.
     */
    template <typename DB>
    void observations(DB *db, Observations &obs);

    /**
     * @brief Add many observations to the database.
     *
     * Use this when adding ~100 or more.
     *
     * @param db An sbsearch database instance.
     *
     * @param obs The observations to add.  `observation_id` must be null and
     *            `terms` must be defined.  `observation_id` will be updated.
     *
     * Raises ObservationError if observation requirements are not met.
     */
    template <typename DB>
    void many_observations(DB *db, Observations &obs);

    /**
     * @brief Add observatory to the database.
     *
     * @param db An sbsearch database instance.
     *
     * @param location The observatory to add.
     *
     */
    template <typename DB>
    void observatory(DB *db, const Observatory &location);

}

#endif // SBSDB_ADD_H_
