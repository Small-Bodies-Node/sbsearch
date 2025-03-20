#ifndef SBSDB_GET_H_
#define SBSDB_GET_H_

#include <cinttypes>
#include <vector>

#include "../observation.h"
#include "../moving_target.h"

using std::vector;

namespace sbsearch::sbsdb::get
{
    /**
     * @brief Get all moving targets in the database.
     *
     * @param db An sbsearch database instance.
     *
     * @return std::vector<MovingTarget>
     */
    template <typename DB>
    vector<MovingTarget> all_moving_targets(DB &db);

    /**
     * @brief Get a moving target by unique moving target ID.
     *
     * @param db An sbsearch database instance.
     *
     * @param moving_target_id The moving target ID.
     *
     * @return MovingTarget
     */
    template <typename DB>
    MovingTarget moving_target(DB &db, int64_t moving_target_id);

    /**
     * @brief Get a moving target by name and small body status.
     *
     * @param db An sbsearch database instance.
     *
     * @param name The target name to search for.
     *
     * @param small_body `true` if the target is a small solar system object.
     *
     * @return MovingTarget The target found in the database, or a new object.
     */
    template <typename DB>
    MovingTarget moving_target(DB &db, const string &name, const bool small_body = true);

    /**
     * @brief Get an observation by ID from the database.
     *
     * @param db An sbsearch database instance.
     *
     * @param observation_id The observation ID.
     *
     * @return Observation
     *
     * Raises ObservationError if the ID is not in the database.
     */
    template <typename DB>
    Observation observation(DB &db, const int64_t observation_id);

    /**
     * @brief Get a set of observations.
     *
     * @param db An sbsearch database instance.
     *
     * @param observation_ids The IDs to get.
     *
     * @return Observations
     *
     * Raises ObservationError if the number of observations returned does not
     * match the length of the vector.
     */
    template <typename DB>
    Observations observations(DB &db, const vector<int64_t> &observation_ids);

}

#endif // SBSDB_GET_H_
