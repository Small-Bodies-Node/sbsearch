#include <algorithm>
#include <cinttypes>
#include <unordered_set>
#include <vector>

#include "add.h"
#include "get.h"
#include "remove.h"
#include "update.h"
#include "verify.h"
#include "postgresql.h"
#include "sbsdb.h"
#include "../exceptions.h"
#include "../moving_target.h"
#include "../observation.h"

using std::endl;
using std::string;
using std::vector;

namespace sbsearch::sbsdb::update
{
    template <typename DB>
    void indexer_options(DB *db, const Indexer::Options &options)
    {
        Logger::info() << "Writing indexer configuration to database." << endl;

        const vector<std::pair<string, int>> parameters{
            {"max_spatial_index_cells", options.max_spatial_index_cells()},
            {"max_spatial_level", options.max_spatial_level()},
            {"min_spatial_level", options.min_spatial_level()},
            {"temporal_resolution", options.temporal_resolution()}};

        for (auto const &parameter : parameters)
            db->template execute(
                "UPDATE configuration SET value = $1 WHERE parameter = $2",
                parameter.second,
                parameter.first);
    }

    template <typename DB>
    void moving_target(DB *db, MovingTarget &target)
    {
        if (!target.moving_target_id())
            throw MovingTargetError("moving_target_id must be defined");

        bool use_transaction = db->template begin();
        try
        {
            // For moving target updates, it is easiest to just remove and add.
            auto old = get::moving_target(db, target.moving_target_id().value());
            remove::moving_target(db, old);
            add::moving_target(db, target);

            if (use_transaction)
                db->template commit();
        }
        catch (const std::exception &err)
        {
            Logger::error() << err.what() << endl;
            if (use_transaction)
                db->template rollback();
        }
    }

    template <typename DB>
    void observations(DB *db, const Observations &obs)
    {
        Logger::info() << "Updating " << obs.size() << " observation"
                       << (obs.size() == 1 ? "" : "s") << "." << endl;

        // observation IDs and terms are required.
        verify::observations(obs, true, true);

        const bool use_transaction = db->template begin();
        try
        {
            for (auto const &observation : obs)
                db->template execute(
                    R"(
                    UPDATE observations
                    SET
                        source=$1,
                        observatory=$2,
                        product_id=$3,
                        mjd_start=$4,
                        mjd_stop=$5,
                        fov=$6,
                        center=$7,
                        terms=$8
                    WHERE observation_id=$9
                )",
                    observation.source(),
                    observation.observatory(),
                    observation.product_id(),
                    observation.mjd_start(),
                    observation.mjd_stop(),
                    observation.fov(),
                    observation.center(),
                    observation.terms(),
                    observation.observation_id().value());

            if (use_transaction)
                db->template commit();
        }
        catch (std::exception &err)
        {
            Logger::error() << err.what() << endl;
            if (use_transaction)
                db->template rollback();
            throw;
        }
    }

    template <typename DB>
    void observations(DB *db, const vector<int64_t> &observation_ids, const vector<vector<string>> &observation_terms)
    {
        Logger::debug() << "Updating index terms for " << observation_ids.size() << " observation"
                        << (observation_ids.size() == 1 ? "" : "s") << "." << endl;

        if (observation_ids.size() != observation_terms.size())
            throw ObservationError("Size of observation_ids does not match size of observation_terms");

        const bool use_transaction = db->template begin();
        try
        {
            for (size_t i = 0; i < observation_ids.size(); i++)
            {
                db->template execute(
                    R"(
                    UPDATE observations
                    SET
                    terms=$1
                    WHERE observation_id=$2
                    )",
                    observation_terms[i],
                    observation_ids[i]);
                cerr << join(observation_terms[i], ",") << endl;
            }

            if (use_transaction)
                db->template commit();
        }
        catch (std::exception &err)
        {
            Logger::error() << err.what() << endl;
            if (use_transaction)
                db->template rollback();
            throw;
        }
    }

    template void indexer_options(Postgresql *, const Indexer::Options &);
    template void moving_target(Postgresql *, MovingTarget &);
    template void observations(Postgresql *, const Observations &);
    template void observations(Postgresql *, const vector<int64_t> &, const vector<vector<string>> &);
}