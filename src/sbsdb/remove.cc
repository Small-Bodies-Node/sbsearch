#include <cinttypes>
#include <vector>

#include "remove.h"
#include "postgresql.h"
#include "../ephemeris.h"
#include "../exceptions.h"
#include "../found.h"
#include "../logging.h"
#include "../moving_target.h"
#include "../observation.h"
#include "../observatory.h"

using std::endl;
using std::string;
using std::vector;

namespace sbsearch::sbsdb::remove
{
    template <typename DB>
    void ephemeris(DB *db, const MovingTarget &target, const double &mjd_start, const double &mjd_stop)
    {
        const int count = db->template get_one<int>(
            "SELECT COUNT(*) FROM ephemerides "
            "WHERE moving_target_id = $1 AND mjd >= $2 AND mjd <= $3",
            target.moving_target_id().value(),
            mjd_start,
            mjd_stop);

        Logger::info() << "Removing " << count
                       << " ephemeris epochs from database for " << target << "." << endl;

        db->template execute(
            "DELETE FROM ephemerides WHERE moving_target_id = $1 AND mjd >= $2 AND mjd <= $3",
            target.moving_target_id(),
            mjd_start,
            mjd_stop);
    }

    template <typename DB>
    void found(DB *db, const Founds &founds)
    {
        Logger::info() << "Removing " << founds.size() << " found entries." << endl;

        const bool use_transaction = db->template begin();
        try
        {
            for (auto const found : founds)
            {
                if (!found.ephemeris.target().moving_target_id())
                    throw MovingTargetError("Cannot remove found items for a moving target without an ID.");

                db->template execute(
                    "DELETE FROM found WHERE moving_target_id=$1 and observation_id=$2",
                    found.ephemeris.target().moving_target_id(),
                    found.observation.observation_id().value());
            }

            if (use_transaction)
                db->template commit();
        }
        catch (const std::exception &err)
        {
            if (use_transaction)
                db->template rollback();
            throw err;
        }
    }

    template <typename DB>
    void moving_target(DB *db, const MovingTarget &target)
    {
        Logger::info() << "Removing " << target << " from the database." << endl;

        if (!target.moving_target_id())
            throw MovingTargetError("moving_target_id is null");

        const int count = db->template get_one<int>(
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
    void observations(DB *db, Observations &observations_)
    {
        Logger::info() << "Removing " << observations_.size() << " observation"
                       << (observations_.size() == 1 ? "" : "s") << "." << endl;

        // Check that observation_id is defined
        verify::observations(observations_, true, false);

        string statement;
        if constexpr (std::is_same_v<DB, Postgresql> == true)
        {
            statement = "WITH deleted AS "
                        "(DELETE FROM observations WHERE observation_id = ANY($1) RETURNING *) "
                        "SELECT COUNT(*) FROM deleted";
        }

        const bool use_transaction = db->template begin();
        try
        {
            int count = db->template get_one<int>(statement, observations_.observation_ids());
            if (count != observations_.size())
                throw ObservationError("only found " + std::to_string(count) + " of " +
                                       std::to_string(observations_.size()) + " observation_ids.");
            if (use_transaction)
                db->template commit();
        }
        catch (std::exception &err)
        {
            if (use_transaction)
                db->template rollback();
            throw;
        }

        // clear the observation IDs
        for (auto it = observations_.begin(); it < observations_.end(); ++it)
            it->observation_id(std::nullopt);
    }

    template <typename DB>
    void observatory(DB *db, const string &name)
    {
        Logger::info() << "Removing observatory with name " << name << endl;
        db->template execute("DELETE FROM observatories WHERE name=$1", name);
    };

    template void ephemeris(Postgresql *, const MovingTarget &, const double &, const double &);
    template void found(Postgresql *, const Founds &);
    template void moving_target(Postgresql *, const MovingTarget &);
    template void observations(Postgresql *, Observations &);
    template void observatory(Postgresql *, const string &);
}
