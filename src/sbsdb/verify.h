#ifndef SBSDB_VERIFY_H_
#define SBSDB_VERIFY_H_

#include "../moving_target.h"
#include "../observation.h"

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
    void moving_target(DB *db, const MovingTarget &target);

    /**
     * @brief Verify that a set of observations are ready for database I/O.
     *
     * @param observations_ The observations to test.
     *
     * @param observation_id_test `true` to test that the observation ID has a
     *                            non-null value, `false` to test that the
     *                            observation ID is null.
     *
     * @param terms_test `true` to test that any terms are present.
     *
     * @throws ObservationError if any observations fail.
     */
    void observations(const Observations &observation_,
                      const bool observation_id_test,
                      const bool terms_test);
}

#endif // SBSDB_VERIFY_H_
