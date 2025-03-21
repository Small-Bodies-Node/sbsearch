#ifndef SBSDB_REMOVE_H_
#define SBSDB_REMOVE_H_

#include <cinttypes>
#include <vector>

// #include "../ephemeris.h"
// #include "../found.h"
// #include "../observation.h"
// #include "../observatory.h"
#include "../moving_target.h"

using std::vector;

namespace sbsearch::sbsdb::remove
{
    /**
     * @brief Remove a moving target from the database.
     *
     * @param db An sbsearch database instance.
     *
     * @param target The moving target (ID) to remove.
     *
     * @return MovingTarget
     */
    template <typename DB>
    void moving_target(DB &db, const MovingTarget &target);

}

#endif // SBSDB_REMOVE_H_
