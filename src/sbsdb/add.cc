#include <algorithm>
#include <cinttypes>
#include <unordered_set>
#include <vector>

#include "./add.h"
#include "./get.h"
#include "./postgresql.h"
#include "./verify.h"
#include "date.h"
#include "ephemeris.h"
#include "exceptions.h"
#include "found.h"
#include "moving_target.h"
#include "observation.h"
#include "observatory.h"

using std::endl;
using std::string;
using std::vector;

namespace sbsearch::sbsdb::add
{
    template <typename DB>
    void ephemeris(DB *db, Ephemeris &ephemeris_)
    {
        Logger::info()
            << "Adding " << ephemeris_.num_vertices()
            << " ephemeris epochs for target " << ephemeris_.target()
            << "." << endl;

        // observation ID and terms are required.
        verify::moving_target(db, ephemeris_.target());

        const bool use_transaction = db->template begin();

        for (const Ephemeris::Datum row : ephemeris_.data())
            db->template execute(
                R"(
                    INSERT INTO ephemerides (
                        moving_target_id, mjd, tmtp,
                        ra, dec, mu, mu_theta,
                        unc_a, unc_b, unc_theta,
                        rh, delta, phase, selong, true_anomaly,
                        sangle, vangle, vmag, mjd_added
                    ) VALUES (
                        $1, $2, $3,
                        $4, $5, $6, $7,
                        $8, $9, $10,
                        $11, $12, $13, $14, $15,
                        $16, $17, $18, $19
                    )
                )",
                ephemeris_.target().moving_target_id(), row.mjd, row.tmtp,
                row.ra, row.dec, row.mu, row.mu_theta,
                row.unc_a, row.unc_b, row.unc_theta,
                row.rh, row.delta, row.phase, row.selong, row.true_anomaly,
                row.sangle, row.vangle, row.vmag, Date::now().mjd());

        if (use_transaction)
            db->template commit();
    }

    template <typename DB>
    void found(DB *db, const Founds &founds)
    {
        Logger::info() << "Adding " << founds.size() << " found observations." << endl;

        const bool use_transaction = db->template begin();

        std::unordered_set<int64_t> checked; // only verify moving_targets once
        try
        {
            for (auto const &found : founds)
            {
                // verify moving_targets, but only once
                int64_t moving_target_id = found.ephemeris.target().moving_target_id().value();
                if (checked.find(moving_target_id) == checked.end())
                    verify::moving_target(db, found.ephemeris.target());
                checked.insert(moving_target_id);

                Ephemeris::Datum eph;
                if (found.ephemeris.num_segments() == 1)
                    eph = found.ephemeris.data(0);
                else
                    eph = found.ephemeris.interpolate(found.observation.mjd_mid()).data(0);

                db->template execute(
                    R"(
                    INSERT INTO found (
                        observation_id, moving_target_id, mjd, tmtp,
                        ra, dec, mu, mu_theta,
                        unc_a, unc_b, unc_theta,
                        rh, delta, phase, selong, true_anomaly,
                        sangle, vangle, vmag, mjd_added
                    ) VALUES (
                        $1, $2, $3, $4,
                        $5, $6, $7, $8,
                        $9, $10, $11,
                        $12, $13, $14, $15, $16,
                        $17, $18, $19, $20
                    )
                )",
                    found.observation.observation_id().value(),
                    found.ephemeris.target().moving_target_id().value(),
                    eph.mjd,
                    eph.tmtp,
                    eph.ra,
                    eph.dec,
                    eph.mu,
                    eph.mu_theta,
                    eph.unc_a,
                    eph.unc_b,
                    eph.unc_theta,
                    eph.rh,
                    eph.delta,
                    eph.phase,
                    eph.selong,
                    eph.true_anomaly,
                    eph.sangle,
                    eph.vangle,
                    eph.vmag,
                    Date::now().mjd());
            }

            if (use_transaction)
            {
                std::cerr << "commitment" << std::endl;
                db->template commit();
            }
        }
        catch (const SBSException &err)
        {
            if (use_transaction)
                db->template rollback();
            throw err;
        }
        catch (const std::exception &err)
        {
            if (use_transaction)
                db->template rollback();
            Logger::error() << err.what() << endl;
            throw err;
        }
    }

    template <typename DB>
    void moving_target(DB *db, MovingTarget &target)
    {
        Logger::info() << "Add moving target " << target << endl;

        // We allow the ID to exist for update::moving_target, but if it does
        // exist, then make sure it is not already in the database
        if (target.moving_target_id())
        {
            MovingTarget test;
            try
            {
                test = get::moving_target(db, target.moving_target_id().value());
            }
            catch (const MovingTargetError &)
            {
                // OK if not found
            }

            if (test.moving_target_id())
                throw MovingTargetError(
                    "moving target id " +
                    std::to_string(target.moving_target_id().value()) +
                    " already exists");
        }

        const bool use_transaction = db->template begin();
        try
        {
            // insert primary designation, getting a new moving_target_id as needed
            int64_t moving_target_id = db->template get_one<int64_t>(
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
                db->template execute(statement,
                                     moving_target_id,
                                     name,
                                     target.small_body());
                Logger::debug() << "Add moving target name " << name
                                << " (ID=" << moving_target_id
                                << "; small body=" << (target.small_body() ? "true" : "false")
                                << "; primary=false)" << std::endl;
            }
            if (use_transaction)
                db->template commit();

            target.moving_target_id(moving_target_id);
        }
        catch (std::exception &err)
        {
            if (use_transaction)
                db->template rollback();
            throw MovingTargetError(err.what());
        }
    }

    template <typename DB>
    void observations(DB *db, Observations &observations_)
    {
        if (observations_.size() > 100)
        {
            many_observations(db, observations_);
            return;
        }

        Logger::info() << "Adding " << observations_.size() << " observation"
                       << (observations_.size() == 1 ? "" : "s") << "." << endl;

        // observation ID must be null, terms are required
        verify::observations(observations_, false, true);

        const bool use_transaction = db->template begin();
        int added = 0;
        try
        {
            for (auto it = observations_.begin(); it < observations_.end(); it++)
            {
                it->mjd_added(Date::now().mjd());
                int64_t observation_id = db->template get_one<int64_t>(
                    R"(
                        INSERT INTO observations
                        (source, observatory, product_id, mjd_start, mjd_stop, fov, 
                         center, terms, meta, mjd_added)
                        VALUES ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10)
                        RETURNING observation_id
                    )",
                    it->source(),
                    it->observatory(),
                    it->product_id(),
                    it->mjd_start(),
                    it->mjd_stop(),
                    it->fov(),
                    it->center(),
                    it->terms(),
                    it->meta(),
                    it->mjd_added());

                it->observation_id(observation_id);
            }
            if (use_transaction)
                db->template commit();
        }
        catch (std::exception &err)
        {
            cerr << err.what() << endl;
            Logger::error() << err.what() << endl;
            if (use_transaction)
                db->template begin();
            throw err;
        }
    }

    template <typename DB>
    void many_observations(DB *db, Observations &observations_)
    {
        Logger::info() << "Adding " << observations_.size() << " observation"
                       << (observations_.size() == 1 ? "" : "s") << "." << endl;

        // observation ID must be null, terms are required
        verify::observations(observations_, false, true);

        const bool use_transaction = db->template begin();
        int added = 0;
        try
        {
            added = db->template insert_many_observations(observations_);
            if (use_transaction)
                db->template commit();
        }
        catch (std::exception &err)
        {
            cerr << err.what() << endl;
            Logger::error() << err.what() << endl;
            if (use_transaction)
                db->template begin();
            throw err;
        }

        Logger::info() << "Added " << added << " observation"
                       << (added == 1 ? "" : "s") << "." << endl;
    }

    template <typename DB>
    void observatory(DB *db, const Observatory &location)
    {
        Logger::info() << "Adding observatory " << location.name << "." << std::endl;

        // do not add anything if this name is already in the database
        if (verify::observatory(db, location.name))
            throw ObservatoryError(location.name + " already exists.");

        db->template execute(
            "INSERT INTO observatories (name, longitude, rho_cos_phi, rho_sin_phi) "
            "VALUES ($1, $2, $3, $4)",
            location.name,
            location.longitude,
            location.rho_cos_phi,
            location.rho_sin_phi);
    }

    template void ephemeris(Postgresql *, Ephemeris &);
    template void found(Postgresql *, const Founds &);
    template void moving_target(Postgresql *, MovingTarget &);
    template void observations(Postgresql *, Observations &);
    template void many_observations(Postgresql *, Observations &);
    template void observatory(Postgresql *, const Observatory &);
}
