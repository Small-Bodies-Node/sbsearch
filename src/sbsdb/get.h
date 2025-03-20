#ifndef SBSDB_GET_H_
#define SBSDB_GET_H_

#include <cinttypes>
#include <vector>

#include "../moving_target.h"

using std::vector;

namespace sbsearch::sbsdb::get
{
    /**
     * @brief Get a moving target by unique moving target ID.
     *
     * @param db An sbsearch database instance.
     * @param id The moving target ID.
     * @return MovingTarget
     */
    template <typename DB>
    MovingTarget moving_target(DB &db, int64_t moving_target_id);

    /**
     * @brief Get a moving target by name and small body status.
     *
     * @param db An sbsearch database instance.
     * @param name The target name to search for.
     * @param small_body `true` if the target is a small solar system object.
     * @return MovingTarget The target found in the database, or a new object.
     */
    template <typename DB>
    MovingTarget moving_target(DB &db, const string &name, const bool small_body = true);

    /// @brief Get all moving targets in the database.
    /// @param db An sbsearch database instance.
    /// @return
    template <typename DB>
    vector<MovingTarget> all_moving_targets(DB &db);
}

#endif // SBSDB_GET_H_
