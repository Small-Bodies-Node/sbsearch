#include "get.h"
#include "postgresql.h"
#include "../exceptions.h"
#include "../moving_target.h"

using std::endl;

namespace sbsearch::sbsdb::verify
{
    template <typename DB>
    void moving_target(DB &db, const MovingTarget &target)
    {
        if (!target.moving_target_id())
            throw MovingTargetError("moving_target_id is not defined");

        // verify that the moving target matches the database
        MovingTarget dbtarget = get::moving_target(
            db, target.moving_target_id().value()); // throws MovingTargetError if not found

        if (target != dbtarget)
            throw MovingTargetError("target does not match database copy");
    }

    template void moving_target(Postgresql &, const MovingTarget &);
}