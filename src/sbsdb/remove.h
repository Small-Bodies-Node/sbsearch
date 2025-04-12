#ifndef SBSDB_REMOVE_H_
#define SBSDB_REMOVE_H_

#include <cinttypes>
#include <vector>

#include "../ephemeris.h"
#include "../found.h"
#include "../moving_target.h"
#include "../observation.h"
#include "../observatory.h"

using std::vector;

namespace sbsearch::sbsdb::remove
{
    /**
     * @brief Remove ephemeris data from the database.
     *
     * @param db An sbsearch database instance.
     *
     * @param target The moving target (ID) ephemeris data to remove.
     *
     * @param mjd_start Remove ephemeris data after this modified Julian date.
     *
     * @param mjd_stop Remove ephemeris data before this modified Julian date.
     */
    template <typename DB>
    void ephemeris(DB *db, const MovingTarget &target, const double &mjd_start = 0, const double &mjd_stop = 100000);

    /**
     * @brief Remove found data from the database.
     *
     * @param db An sbsearch database instance.
     *
     * @param founds The found data to remove.
     */
    template <typename DB>
    void found(DB *db, const Founds &founds);

    /**
     * @brief Remove a moving target from the database.
     *
     * @param db An sbsearch database instance.
     *
     * @param target The moving target (ID) to remove.
     */
    template <typename DB>
    void moving_target(DB *db, const MovingTarget &target);

    /**
     * @brief Remove observations from the database.
     *
     * @param db An sbsearch database instance.
     *
     * @param observations_ The observations to remove.
     *
     * @returns Observations A copy of the observations with the observation IDs
     *          removed.
     */
    template <typename DB>
    void observations(DB *db, Observations &observations_);

    /**
     * @brief Remove an observatory from the database.
     *
     * @param db An sbsearch database instance.
     *
     * @param name The name of the observatory to remove.
     */
    template <typename DB>
    void observatory(DB *db, const string &name);

}

#endif // SBSDB_REMOVE_H_
