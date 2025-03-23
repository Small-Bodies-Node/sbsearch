#ifndef SBSDB_VERIFY_H_
#define SBSDB_VERIFY_H_

#include "../moving_target.h"

namespace sbsearch::sbsdb::verify
{
    /**
     * @brief Verify that the moving target exists and metadata is up to date.
     *
     * @param db An sbsearch database.
     *
     * @param target The target to check.
     *
     * @throws MovingTargetError If the target does not match the database.
     */
    template <typename DB>
    void moving_target(DB &db, const MovingTarget &target);
}

#endif // SBSDB_VERIFY_H_
