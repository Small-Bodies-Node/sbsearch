#ifndef SBSDB_GET_H_
#define SBSDB_GET_H_

#include <cinttypes>
#include <optional>
#include <utility>
#include <vector>

#include "../ephemeris.h"
#include "../found.h"
#include "../indexer.h"
#include "../observation.h"
#include "../observatory.h"
#include "../moving_target.h"

using std::vector;

namespace sbsearch::sbsdb::get
{
    /**
     * @brief Get all moving targets in the database.
     *
     * @param db An sbsearch database instance.
     *
     * @return std::vector<MovingTarget>
     */
    template <typename DB>
    vector<MovingTarget> all_moving_targets(DB &db);

    /**
     * @brief Get all observatories in the database.
     *
     * @param db An sbsearch database instance.
     *
     * @return Observatories
     */
    template <typename DB>
    Observatories all_observatories(DB &db);

    /**
     * @brief Get a moving target's ephemeris from the database.
     *
     * @param db An sbsearch database instance.
     *
     * @param target The moving target.
     *
     * @return Ephemeris
     */
    template <typename DB>
    Ephemeris ephemeris(DB &db,
                        const MovingTarget &target,
                        double mjd_start = 0,
                        double mjd_stop = 100000);

    /**
     * @brief Get the date range for a target's ephemeris.
     *
     * @param db An sbsearch database instance.
     *
     * @param target The moving target to consider.
     *
     * @return std::pair<optional<double>, optional<double>> The date range as a
     *                                                       pair of modified
     *                                                       Julian dates.
     */
    template <typename DB>
    std::pair<optional<double>, optional<double>>
    ephemeris_date_range(DB &db, const MovingTarget &target);

    /**
     * @brief Get found object data.
     *
     * @param db An sbsearch database instance.
     *
     * @param observation Find found observations for this target.
     *
     * @return Founds
     */
    template <typename DB>
    Founds found(DB &db, const MovingTarget &target);

    /**
     * @brief Get found object data.
     *
     * @param db An sbsearch database instance.
     *
     * @param observation Find found targets for this observation.
     *
     * @return Founds
     */
    template <typename DB>
    Founds found(DB &db, const Observation &observation);

    /**
     * @brief Get the indexer options from the database.
     *
     * @param db An sbsearch database instance.
     *
     * @return Indexer::Options
     */
    template <typename DB>
    Indexer::Options indexer_options(DB &db);

    /**
     * @brief Get a moving target by unique moving target ID.
     *
     * @param db An sbsearch database instance.
     *
     * @param moving_target_id The moving target ID.
     *
     * @return MovingTarget
     */
    template <typename DB>
    MovingTarget moving_target(DB &db, int64_t moving_target_id);

    /**
     * @brief Get a moving target by name and small body status.
     *
     * @param db An sbsearch database instance.
     *
     * @param name The target name to search for.
     *
     * @param small_body `true` if the target is a small solar system object.
     *
     * @return MovingTarget The target found in the database, or a new object.
     */
    template <typename DB>
    MovingTarget moving_target(DB &db, const string &name, const bool small_body = true);

    /**
     * @brief Get the date range for the observations table.
     *
     * @param db An sbsearch database instance.
     *
     * @param source Optionally limit the range to this data source.
     *
     * @return std::pair<optional<double>, optional<double>> The date range as a
     *                                                       pair of modified
     *                                                       Julian dates.
     */
    template <typename DB>
    std::pair<optional<double>, optional<double>>
    observations_date_range(DB &db, const optional<string> &source = {});

    /**
     * @brief Get a set of observations.
     *
     * @param db An sbsearch database instance.
     *
     * @param observation_ids The IDs to get.
     *
     * @return Observations
     *
     * Raises ObservationError if the number of observations returned does not
     * match the length of the vector.
     */
    template <typename DB>
    Observations observations(DB &db, const vector<optional<int64_t>> &observation_ids);

    /**
     * @brief Get an observatory by name from the database.
     *
     * @param db An sbsearch database instance.
     *
     * @param name The observatory name.
     *
     * @return Observatory
     *
     * Raises ObservatoryError if the name is not in the database.
     */
    template <typename DB>
    Observatory observatory(DB &db, const string &name);

    /**
     * @brief Get a list of all sources.
     *
     * @param db An sbsearch database instance.
     *
     * @param name The observatory name.
     *
     * @return Observatory
     *
     * Raises ObservatoryError if the name is not in the database.
     */
    template <typename DB>
    vector<string> sources(DB &db);
}

#endif // SBSDB_GET_H_
