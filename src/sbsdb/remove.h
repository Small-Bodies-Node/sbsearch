#ifndef SBSDB_REMOVE_H_
#define SBSDB_REMOVE_H_

#include <cinttypes>
#include <vector>

// #include "../ephemeris.h"
// #include "../found.h"
#include "../moving_target.h"
// #include "../observation.h"
#include "../observatory.h"

using std::vector;

namespace sbsearch::sbsdb::remove
{
    /**
     * @brief Remove a moving target from the database.
     *
     * @param db An sbsearch database instance.
     *
     * @param target The moving target (ID) to remove.
     */
    template <typename DB>
    void moving_target(DB &db, const MovingTarget &target);

    /**
     * @brief Remove an observatory from the database.
     *
     * @param db An sbsearch database instance.
     *
     * @param name The name of the observatory to remove.
     */
    template <typename DB>
    void observatory(DB &db, const string &name);

}

#endif // SBSDB_REMOVE_H_
