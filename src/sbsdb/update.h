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
     * @brief Write indexer options to the database.
     *
     * @param db An sbsearch database instance.
     *
     * @param options The options to write.
     */
    template <typename DB>
    void indexer_options(DB *db, const Indexer::Options &options);

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
    void moving_target(DB *db, MovingTarget &target);

    /**
     * @brief Update observations in the database.
     *
     * @param db An sbsearch database instance.
     *
     * @param obs The observations to update.  The database copies are
     *            identified by `product_id`. `observation_id` may be updated.
     *
     * Raises ObservationError if requirements are not met.
     */
    template <typename DB>
    void observations(DB *db, Observations &obs);

    /**
     * @brief Specialized updater for re-indexing observations.
     *
     * Best for updating many observations (>>100).
     *
     * @param db An sbsearch database instance.
     *
     * @param observation_ids The observations to update.
     *
     * @param observation_terms The terms for each observation.
     */
    template <typename DB>
    void observations(DB *db, const vector<int64_t> &observation_ids, const vector<vector<string>> &observation_terms);
}

#endif // SBSDB_UPDATE_H_
