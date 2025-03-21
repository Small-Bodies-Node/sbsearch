#include <algorithm>
#include <cinttypes>
#include <unordered_set>
#include <vector>

#include "add.h"
#include "get.h"
#include "remove.h"
#include "update.h"
#include "postgresql.h"
#include "sbsdb.h"
// #include "../ephemeris.h"
// #include "../exceptions.h"
// #include "../found.h"
#include "../moving_target.h"
// #include "../observation.h"
// #include "../observatory.h"

using std::endl;
using std::string;
using std::vector;

namespace sbsearch::sbsdb::update
{
    template <typename DB>
    void moving_target(DB &db, MovingTarget &target)
    {
        if (!target.moving_target_id())
            throw MovingTargetError("moving_target_id must be defined");

        bool use_transaction = !(db.template in_transaction());
        if (use_transaction)
            db.template begin();

        try
        {
            // For moving target updates, it is easiest to just remove and add.
            auto old = get::moving_target(db, target.moving_target_id().value());
            remove::moving_target(db, old);
            add::moving_target(db, target);

            if (use_transaction)
                db.template commit();
        }
        catch (const std::exception &err)
        {
            Logger::error() << err.what() << endl;
            if (use_transaction)
                db.template rollback();
        }
    }

    template void moving_target(Postgresql &, MovingTarget &);
}