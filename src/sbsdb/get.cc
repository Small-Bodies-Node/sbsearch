#include <algorithm>
#include <cinttypes>
#include <optional>
#include <unordered_set>
#include <utility>
#include <vector>

#include "get.h"
#include "postgresql.h"
#include "sbsdb.h"
#include "../ephemeris.h"
#include "../exceptions.h"
#include "../found.h"
#include "../indexer.h"
#include "../moving_target.h"
#include "../observation.h"
#include "../observatory.h"

using std::endl;
using std::string;
using std::vector;

namespace sbsearch::sbsdb::get
{
    template <typename DB>
    vector<MovingTarget> all_moving_targets(DB *db)
    {
        auto rows = db->template get_many<MovingTarget::DBModel>(
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
                [&target](const MovingTarget::DBModel &row)
                { return row.moving_target_id == target.moving_target_id(); });

            std::for_each(
                start,
                end,
                [&target](const MovingTarget::DBModel &row)
                { target.add_name(row.name, row.primary_id); });

            result.push_back(target);
            start = end;
        } while (start < rows.end());

        return result;
    }

    template <typename DB>
    Observatories all_observatories(DB *db)
    {
        auto obs_vector = db->template get_many<Observatory>("SELECT * FROM observatories");
        Observatories observatories;
        for (auto &observatory : obs_vector)
            observatory.insert_into(observatories);
        return observatories;
    }

    template <typename DB>
    vector<std::pair<int64_t, string>> all_observations_fov(DB *db, const int limit, const int64_t offset)
    {
        return db->template get_many<int64_t, string>(
            "SELECT observation_id,fov FROM observations LIMIT $1 OFFSET $2",
            limit, offset);
    }

    template <typename DB>
    Ephemeris ephemeris(DB *db, const MovingTarget &target, const double mjd_start, const double mjd_stop)
    {
        if (!target.moving_target_id())
            throw MovingTargetError("Cannot get ephemeris for moving target with an undefined ID.");

        const int count = db->template get_one<int>(
            "SELECT COUNT(*) FROM ephemerides "
            "WHERE moving_target_id=$1 AND mjd >= $2 and mjd <= $3",
            target.moving_target_id().value(),
            mjd_start,
            mjd_stop);

        Logger::debug() << "Reading " << count
                        << " ephemeris epochs from database for " << target.designation()
                        << " (moving_target_id=" << target.moving_target_id().value() << ")"
                        << endl;

        Ephemeris::Data data(
            db->template get_many<Ephemeris::Datum>(
                "SELECT * FROM ephemerides "
                "WHERE moving_target_id = $1 AND mjd >= $2 and mjd <= $3",
                target.moving_target_id().value(),
                mjd_start,
                mjd_stop));

        return {target, data};
    }

    template <typename DB>
    std::pair<optional<double>, optional<double>>
    ephemeris_date_range(DB *db, const MovingTarget &target)
    {
        if (!target.moving_target_id())
            throw MovingTargetError("moving_target_id is null");

        auto mjd_start = db->template get_one<optional<double>>(
            "SELECT MIN(mjd) FROM ephemerides WHERE moving_target_id=$1",
            target.moving_target_id().value());
        auto mjd_stop = db->template get_one<optional<double>>(
            "SELECT MAX(mjd) FROM ephemerides WHERE moving_target_id=$1",
            target.moving_target_id().value());

        return {mjd_start, mjd_stop};
    }

    template <typename DB>
    Founds found(DB *db, const MovingTarget &target)
    {
        if (!target.moving_target_id())
            throw MovingTargetError("moving_target_id is null");

        auto rows = db->template get_many<Found::DBModel>(
            "SELECT * FROM found WHERE moving_target_id=$1",
            target.moving_target_id().value());

        Founds result;
        for (const Found::DBModel &row : rows)
        {
            Observation obs = observations(db, {row.observation_id})[0];
            Ephemeris ephemeris(
                target, {{
                            row.mjd,
                            row.tmtp,
                            row.ra,
                            row.dec,
                            row.unc_a,
                            row.unc_b,
                            row.unc_theta,
                            row.rh,
                            row.delta,
                            row.phase,
                            row.selong,
                            row.true_anomaly,
                            row.sangle,
                            row.vangle,
                            row.vmag,
                        }});
            result.append(Found(obs, ephemeris));
        }
        return result;
    }

    template <typename DB>
    Founds found(DB *db, const Observation &obs)
    {
        auto rows = db->template get_many<Found::DBModel>(
            "SELECT * FROM found WHERE observation_id=$1",
            obs.observation_id());

        Founds result;
        for (const Found::DBModel &row : rows)
        {
            MovingTarget target = moving_target(db, row.moving_target_id);
            Ephemeris ephemeris(
                target, {{
                            row.mjd,
                            row.tmtp,
                            row.ra,
                            row.dec,
                            row.unc_a,
                            row.unc_b,
                            row.unc_theta,
                            row.rh,
                            row.delta,
                            row.phase,
                            row.selong,
                            row.true_anomaly,
                            row.sangle,
                            row.vangle,
                            row.vmag,
                        }});
            result.append(Found(obs, ephemeris));
        }
        return result;
    }

    template <typename DB>
    Indexer::Options indexer_options(DB *db)
    {
        Indexer::Options options;

        options.max_spatial_index_cells(
            db->template get_one<int>(
                "SELECT value FROM configuration WHERE parameter='max_spatial_index_cells'"));
        options.max_spatial_level(
            db->template get_one<int>(
                "SELECT value FROM configuration WHERE parameter='max_spatial_level'"));
        options.min_spatial_level(
            db->template get_one<int>(
                "SELECT value FROM configuration WHERE parameter='min_spatial_level'"));
        options.temporal_resolution(
            db->template get_one<int>(
                "SELECT value FROM configuration WHERE parameter='temporal_resolution'"));

        return options;
    }

    template <typename DB>
    MovingTarget moving_target(DB *db, const int64_t moving_target_id)
    {
        MovingTarget result;
        result.moving_target_id(moving_target_id);

        // one target per name
        auto rows = db->template get_many<MovingTarget::DBModel>(
            "SELECT * FROM moving_targets WHERE moving_target_id = $1", moving_target_id);

        if (rows.size() == 0)
            throw MovingTargetError("moving target id " +
                                    std::to_string(moving_target_id) +
                                    " not found");

        // Package the results into a single target.
        result.small_body(rows[0].small_body);
        for (auto const &row : rows)
            result.add_name(row.name, row.primary_id);

        return result;
    }

    template <typename DB>
    MovingTarget moving_target(DB *db, const string &name, const bool small_body)
    {
        MovingTarget result(name, small_body);

        auto rows = db->template get_many<MovingTarget::DBModel>(
            "SELECT * FROM moving_targets WHERE moving_target_id = "
            "(SELECT moving_target_id FROM moving_targets WHERE name=$1 AND small_body=$2)",
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
    std::pair<optional<double>, optional<double>>
    observations_date_range(DB *db, const optional<string> &source)
    {
        optional<double> mjd_start;
        optional<double> mjd_stop;

        if (!source)
        {
            mjd_start = db->template get_one<optional<double>>(
                "SELECT MIN(mjd_start) FROM observations");
            mjd_stop = db->template get_one<optional<double>>(
                "SELECT MAX(mjd_stop) FROM observations");
        }
        else
        {
            mjd_start = db->template get_one<optional<double>>(
                "SELECT MIN(mjd_start) FROM observations WHERE source=$1", source);
            mjd_stop = db->template get_one<optional<double>>(
                "SELECT MAX(mjd_stop) FROM observations WHERE source=$1", source);
        }

        return {mjd_start, mjd_stop};
    }

    template <typename DB>
    Observations observations(DB *db, const vector<optional<int64_t>> &observation_ids)
    {
        const int count = std::count_if(observation_ids.begin(), observation_ids.end(),
                                        [](auto const &obs)
                                        { return !obs.has_value(); });
        if (count > 0)
            throw ObservationError(std::to_string(count) + " null observation IDs");

        string statement = "SELECT * FROM observations WHERE ";
        if constexpr (std::is_same_v<DB, Postgresql> == true)
            statement += "observation_id = ANY($1)";
        else
            statement += "observation_id IN $1";

        vector<Observation> results = db->template get_many<Observation>(statement, observation_ids);

        if (results.size() != observation_ids.size())
            throw ObservationError(
                "Only found " + std::to_string(results.size()) + " of " +
                std::to_string(observation_ids.size()) + " observations.");

        return {results};
    }

    template <typename DB>
    Observatory observatory(DB *db, const string &name)
    {
        const int count = db->template get_one<int>(
            "SELECT COUNT(*) FROM observatories WHERE name=$1",
            name);
        if (count == 0)
            throw ObservatoryError(name + " not found");

        try
        {
            return db->template get_one<Observatory>(
                "SELECT * FROM observatories WHERE name=$1",
                name);
        }
        catch (std::exception &err)
        {
            throw ObservatoryError(err.what());
        }
    }

    template <typename DB>
    vector<string> sources(DB *db)
    {
        return db->template get_many<string>("SELECT DISTINCT(source) FROM observations");
    }

    template vector<MovingTarget> all_moving_targets(Postgresql *);
    template Observatories all_observatories(Postgresql *);
    template vector<std::pair<int64_t, string>> all_observations_fov(Postgresql *, const int, const int64_t);
    template Ephemeris ephemeris(Postgresql *, const MovingTarget &, const double, const double);
    template std::pair<optional<double>, optional<double>> ephemeris_date_range(Postgresql *, const MovingTarget &);
    template Founds found(Postgresql *, const MovingTarget &);
    template Founds found(Postgresql *, const Observation &);
    template Indexer::Options indexer_options(Postgresql *);
    template MovingTarget moving_target(Postgresql *, const int64_t);
    template MovingTarget moving_target(Postgresql *, const string &, const bool);
    template Observations observations(Postgresql *, const vector<optional<int64_t>> &);
    template std::pair<optional<double>, optional<double>> observations_date_range(Postgresql *, const optional<string> &);
    template Observatory observatory(Postgresql *, const string &);
    template vector<string> sources(Postgresql *);
}