#include <cinttypes>
#include <vector>

#include "remove.h"
#include "postgresql.h"
// #include "sbsdb.h"
// #include "../ephemeris.h"
// #include "../exceptions.h"
// #include "../found.h"
#include "../logging.h"
#include "../moving_target.h"
// #include "../observation.h"
// #include "../observatory.h"

using std::endl;
using std::string;
using std::vector;

namespace sbsearch::sbsdb::remove
{
    template <typename DB>
    void moving_target(DB &db, const MovingTarget &target)
    {
        Logger::info() << "Removing " << target << " from the database." << endl;

        if (!target.moving_target_id())
            throw MovingTargetError("moving_target_id is null");

        const int count = db.template get_one<int>(
            R"(
                    WITH deleted AS
                    (DELETE FROM moving_targets WHERE moving_target_id=$1
                    RETURNING *)
                    SELECT COUNT(*) FROM deleted;
                )",
            target.moving_target_id());

        if (count == 0)
            Logger::warning() << "moving_target_id="
                              << target.moving_target_id().value()
                              << " not found." << endl;

        Logger::info() << count << " moving target names deleted." << endl;
    }

    template void moving_target(Postgresql &, const MovingTarget &);
}
