#include <algorithm>
#include <cinttypes>
#include <unordered_set>
#include <vector>

#include "add.h"
#include "get.h"
#include "postgresql.h"
#include "sbsdb.h"
#include "../ephemeris.h"
#include "../exceptions.h"
#include "../found.h"
#include "../moving_target.h"
#include "../observation.h"
#include "../observatory.h"

using std::endl;
using std::string;
using std::vector;

namespace sbsearch::sbsdb::add
{
    template <typename DB>
    void moving_target(DB &db, MovingTarget &target)
    {
        Logger::info() << "Add moving target " << target << endl;

        // We allow the ID to exist for update::moving_target, but if it does
        // exist, then make sure it is not already in the database
        if (target.moving_target_id())
        {
            const int count = db.template get_one<int>(
                "SELECT COUNT(*) FROM moving_targets WHERE moving_target_id=$1",
                target.moving_target_id());

            if (count != 0)
                throw MovingTargetError("moving target id " +
                                        std::to_string(target.moving_target_id().value()) +
                                        " already exists");
        }

        bool use_transaction = !(db.template in_transaction());
        try
        {
            if (use_transaction)
                db.template begin();

            // insert primary designation, getting a new moving_target_id as needed
            int64_t moving_target_id = db.template get_one<int64_t>(
                R"(
                INSERT INTO moving_targets
                (moving_target_id, name, small_body, primary_id)
                SELECT
                    COALESCE($1, MAX(moving_target_id) + 1, 1) AS moving_target_id,
                    $2 AS name,
                    $3 AS small_body,
                    TRUE AS primary_id
                FROM moving_targets
                RETURNING moving_target_id;
            )",
                target.moving_target_id(),
                target.designation(),
                target.small_body());

            Logger::debug() << "Add moving target name " << target.designation()
                            << " (ID=" << moving_target_id
                            << "; small body=" << (target.small_body() ? "true" : "false")
                            << "; primary=true)" << std::endl;

            // insert alternate names
            const string statement = "INSERT INTO moving_targets "
                                     "(moving_target_id, name, small_body, primary_id) "
                                     "VALUES ($1, $2, $3, FALSE)";

            for (const string &name : target.alternate_names())
            {
                db.template execute(statement,
                                    moving_target_id,
                                    name,
                                    target.small_body());
                Logger::debug() << "Add moving target name " << name
                                << " (ID=" << moving_target_id
                                << "; small body=" << (target.small_body() ? "true" : "false")
                                << "; primary=false)" << std::endl;
            }
            if (use_transaction)
                db.template commit();

            target.moving_target_id(moving_target_id);
        }
        catch (std::exception &err)
        {
            if (use_transaction)
                db.template rollback();
            throw MovingTargetError(err.what());
        }
    }

    template <typename DB>
    Observations observations(DB &db, const Observations &obs)
    {
        Logger::info() << "Adding " << obs.size() << " observation"
                       << (obs.size() == 1 ? "" : "s") << "." << endl;

        int count = std::count_if(obs.begin(), obs.end(),
                                  [](auto const &o)
                                  { return o.observation_id().has_value(); });
        if (count != 0)
            throw ObservationError(std::to_string(count) +
                                   " observations have non-null observation_id");

        count = std::count_if(obs.begin(), obs.end(),
                              [](auto const &o)
                              { return o.terms().size() == 0; });
        if (count != 0)
            throw ObservationError(std::to_string(count) +
                                   " observations missing terms");

        int added = 0;
        Observations result(obs);
        db.template execute("BEGIN");
        try
        {
            for (vector<Observation>::iterator it = result.begin(); it < result.end(); it++)
            {
                int64_t observation_id = db.template get_one<int64_t>(
                    R"(
                        INSERT INTO observations
                        (source, observatory, product_id, mjd_start, mjd_stop, fov, terms)
                        VALUES ($1, $2, $3, $4, $5, $6, $7)
                        RETURNING observation_id
                    )",
                    it->source(),
                    it->observatory(),
                    it->product_id(),
                    it->mjd_start(),
                    it->mjd_stop(),
                    it->fov(),
                    it->terms());

                it->observation_id(observation_id);
            }
            db.template execute("COMMIT");
        }
        catch (std::exception &err)
        {
            Logger::error() << err.what() << endl;
            db.template execute("ROLLBACK");
            throw err;
        }

        return result;
    }

    template <typename DB>
    void observatory(DB &db, const Observatory &location)
    {
        Logger::info() << "Adding observatory " << location.name << "." << std::endl;

        // do not add anything if this name is already in the database
        try
        {
            get::observatory(db, location.name);
        }
        catch (const ObservatoryError &e)
        {
            db.template execute(
                "INSERT INTO observatories (name, longitude, rho_cos_phi, rho_sin_phi) "
                "VALUES ($1, $2, $3, $4)",
                location.name,
                location.longitude,
                location.rho_cos_phi,
                location.rho_sin_phi);
            return;
        }

        throw ObservatoryError(location.name + " already exists.");
    }

    template void moving_target(Postgresql &, MovingTarget &);
    template Observations observations(Postgresql &, const Observations &);
    template void observatory(Postgresql &, const Observatory &);
}
