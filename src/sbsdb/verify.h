#ifndef SBSDB_VERIFY_H_
#define SBSDB_VERIFY_H_

#include <string_view>

#include "moving_target.h"
#include "observation.h"

using std::string_view;

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
     * @param terms_test `true` to test that index terms are present.
     *
     * @throws ObservationError if any observations fail.
     */
    void observations(const Observations &observation_,
                      const bool observation_id_test,
                      const bool terms_test);

    /**
     * @brief Verify that an observatory exists in the database.
     *
     * @param db An sbsearch database.
     *
     * @param name The observatory name to test.
     */
    template <typename DB>
    bool observatory(DB *db, string_view name);
}

#endif // SBSDB_VERIFY_H_
