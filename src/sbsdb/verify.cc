#include "get.h"
#include "postgresql.h"
#include "../exceptions.h"
#include "../moving_target.h"
#include "../observation.h"

using std::endl;

namespace sbsearch::sbsdb::verify
{
    template <typename DB>
    void moving_target(DB *db, const MovingTarget &target)
    {
        if (!target.moving_target_id())
            throw MovingTargetError("moving_target_id is not defined");

        // verify that the moving target matches the database
        MovingTarget dbtarget = get::moving_target(
            db, target.moving_target_id().value()); // throws MovingTargetError if not found

        if (target != dbtarget)
            throw MovingTargetError("target does not match database copy");
    }

    void observations(const Observations &observations_,
                      const bool observation_id_test,
                      const bool terms_test)
    {
        int count = std::count_if(
            observations_.begin(), observations_.end(),
            [&](auto const &observation)
            { return observation.observation_id().has_value() != observation_id_test; });

        if (count != 0)
            throw ObservationError(std::to_string(count) +
                                   " observation ID(s) are " +
                                   (observation_id_test ? "" : "not") +
                                   " null");

        if (terms_test)
        {
            count = std::count_if(observations_.begin(), observations_.end(),
                                  [](auto const &observation)
                                  { return observation.terms().size() == 0; });
            if (count != 0)
                throw ObservationError(std::to_string(count) +
                                       " observations missing terms");
        }
    }

    template void moving_target(Postgresql *, const MovingTarget &);
}