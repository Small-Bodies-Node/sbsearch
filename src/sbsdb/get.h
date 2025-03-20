#ifndef SBSDB_GET_H_
#define SBSDB_GET_H_

#include <cinttypes>
#include <vector>

#include "../moving_target.h"

using std::vector;

namespace sbsearch::sbsdb
{
    // Get a moving target by database ID
    template <typename DB>
    MovingTarget get_moving_target(DB &db, int64_t id);

    /**
     * @brief Get a moving target by name and small body status.
     *
     * @tparam DB An sbsearch::sbsdb class.
     * @param db An sbsearch database instance.
     * @param name The target name to search for.
     * @param small_body `true` if the target is a small solar system object.
     * @return MovingTarget The target found in the database, or a new object.
     */
    template <typename DB>
    MovingTarget get_moving_target(DB &db, const string &name, const bool small_body = true);

    /// @brief
    /// @tparam DB
    /// @param db
    /// @return
    template <typename DB>
    vector<MovingTarget> get_all_moving_targets(DB &db);
}

#endif // SBSDB_GET_H_
