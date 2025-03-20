
#include <algorithm>
#include <cinttypes>
#include <unordered_set>
#include <vector>

#include "get.h"
#include "postgresql.h"
#include "sbsdb.h"
#include "../exceptions.h"
#include "../moving_target.h"
#include "../observation.h"

using std::vector;

namespace sbsearch::sbsdb::get
{
    template <typename DB>
    vector<MovingTarget> all_moving_targets(DB &db)
    {

        auto rows = db.template get_many<MovingTargetsModel>(
            "SELECT * FROM moving_targets ORDER BY moving_target_id");

        vector<MovingTarget> result;

        // Iterate over all rows, identifying targets by moving_target_id
        auto start = rows.begin();
        do
        {
            MovingTarget target;
            target.moving_target_id(start->moving_target_id);
            target.small_body(start->small_body);

            auto end = std::find_if_not(
                start,
                rows.end(),
                [&target](const MovingTargetsModel &row)
                { return row.moving_target_id == target.moving_target_id(); });

            std::for_each(
                start,
                end,
                [&target](const MovingTargetsModel &row)
                { target.add_name(row.name, row.primary_id); });

            result.push_back(target);
            start = end;
        } while (start < rows.end());

        return result;
    }

    template <typename DB>
    MovingTarget moving_target(DB &db, int64_t moving_target_id)
    {
        MovingTarget result;
        result.moving_target_id(moving_target_id);

        // one target per name
        auto rows = db.template get_many<MovingTargetsModel>(
            "SELECT * FROM moving_targets WHERE moving_target_id = $1", moving_target_id);

        // Package the results into a single target.
        result.small_body(rows[0].small_body);
        for (auto const &row : rows)
            result.add_name(row.name, row.primary_id);

        return result;
    }

    template <typename DB>
    MovingTarget moving_target(DB &db, const string &name, const bool small_body)
    {
        MovingTarget result(name, small_body);

        auto rows = db.template get_many<MovingTargetsModel>(
            "SELECT * FROM moving_targets WHERE moving_target_id = "
            "(SELECT moving_target_id FROM moving_targets WHERE name = $1 AND small_body = $2)",
            name, small_body);

        if (rows.size() == 0)
            // this name-small body / name combo is not in the database, return a new object
            return MovingTarget(name, small_body);

        // Otherwise, package the results into a single target.
        result.moving_target_id(rows[0].moving_target_id);
        for (auto const &row : rows)
            result.add_name(row.name, row.primary_id);

        return result;
    }

    template <typename DB>
    Observation observation(DB &db, const int64_t observation_id)
    {
        Observation obs = db.template get_one<Observation>(
            "SELECT * FROM observations WHERE observation_id=$1",
            observation_id);
    }

    template vector<MovingTarget> all_moving_targets(Postgresql &);
    template MovingTarget moving_target(Postgresql &, int64_t);
    template MovingTarget moving_target(Postgresql &, const string &, const bool);
}