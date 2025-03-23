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
    void indexer_options(DB &db, const Indexer::Options &options)
    {
        Logger::info() << "Writing indexer configuration to database." << endl;

        vector<pair<string, int>> parameters{
            {"max_spatial_index_cells", options.max_spatial_index_cells()},
            {"max_spatial_level", options.max_spatial_level()},
            {"min_spatial_level", options.min_spatial_level()},
            {"temporal_resolution", options.temporal_resolution()}};

        for (auto const &parameter : parameters)
            db.template execute(
                "UPDATE configuration SET value = $1 WHERE parameter = $2",
                parameter.first,
                parameter.second);
    }

    template <typename DB>
    void moving_target(DB &db, MovingTarget &target)
    {
        if (!target.moving_target_id())
            throw MovingTargetError("moving_target_id must be defined");

        bool use_transaction = db.template begin();
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

    template <typename DB>
    void observations(DB &db, const Observations &obs)
    {
        Logger::info() << "Updating " << obs.size() << " observation"
                       << (obs.size() == 1 ? "" : "s") << "." << endl;

        // observation IDs and terms are required.
        verify::observations(obs, true, true);

        const bool use_transaction = db.template begin();
        for (auto const &observation : obs)
            db.template execute(
                R"(
                    UPDATE observations
                    SET
                        source=$1,
                        observatory=$2,
                        product_id=$3,
                        mjd_start=$4,
                        mjd_stop=$5,
                        fov=$6,
                        terms=$7
                    WHERE observation_id=$8
                )",
                observation.source(),
                observation.observatory(),
                observation.product_id(),
                observation.mjd_start(),
                observation.mjd_stop(),
                observation.fov(),
                observation.terms(),
                observation.observation_id().value());
    }

    template void moving_target(Postgresql &, MovingTarget &);
    template void observations(Postgresql &, const Observations &);
}