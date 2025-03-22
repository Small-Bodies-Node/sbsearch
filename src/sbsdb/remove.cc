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
    void ephemeris(DB &db, const MovingTarget &target, const double &mjd_start, const double &mjd_stop)
    {
        const int count = db.template get_one<int>(
            "SELECT COUNT(*) FROM ephemerides "
            "WHERE moving_target_id = $1 AND mjd >= $2 AND mjd <= $3",
            target.moving_target_id().value(),
            mjd_start,
            mjd_stop);

        Logger::info() << "Removing " << count
                       << " ephemeris epochs from database for " << target << "." << endl;

        db.template execute(
            "DELETE FROM ephemerides WHERE moving_target_id = $1 AND mjd >= $2 AND mjd <= $3",
            target.moving_target_id(),
            mjd_start,
            mjd_stop);
    }

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

    template <typename DB>
    void observatory(DB &db, const string &name)
    {
        Logger::info() << "Removing observatory with name " << name << endl;
        db.template execute("DELETE FROM observatories WHERE name=$1", name);
    };

    template void ephemeris(Postgresql &, const MovingTarget &, const double &, const double &);
    template void moving_target(Postgresql &, const MovingTarget &);
    template void observatory(Postgresql &, const string &);
}
