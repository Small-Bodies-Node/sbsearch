#ifndef SBSDB_H_
#define SBSDB_H_

#include <string>
#include <utility>
#include <s2/s2point.h>
#include <s2/s2cap.h>
#include <s2/s2metrics.h>
#include <s2/s2polyline.h>
#include <s2/s2polygon.h>
#include <s2/s2region_term_indexer.h>

#include "ephemeris.h"
#include "found.h"
#include "indexer.h"
#include "moving_target.h"
#include "observation.h"
#include "observatory.h"

#define SBSEARCH_DATABASE_VERSION "3.0"

using std::string;

namespace sbsearch
{
    class SBSearchDatabase
    {
    public:
        // close database connection
        virtual void close() = 0;

        // initialize database, or add any missing tables, indices, etc.
        // must be idempotent
        virtual void setup_tables() = 0;

        // execute a sql statement
        virtual void execute_sql(const char *statement) const = 0;

        // get single value results from a SQL statement
        virtual double *get_double(const char *statement) const = 0;
        virtual int *get_int(const char *statement) const = 0;
        virtual int64 *get_int64(const char *statement) const = 0;
        virtual string *get_string(const char *statement) const = 0;

        // drop/create observations indices, e.g., before inserting many new observations
        // restore indices with setup_tables()
        virtual void drop_observations_indices() = 0;
        virtual void create_observations_indices() = 0;

        // Read indexer options from the database.
        Indexer::Options indexer_options();
        // Write indexer options to the database.
        virtual void indexer_options(Indexer::Options options) = 0;

        // get date range, optionally for a single source
        virtual std::pair<double *, double *> observation_date_range(const string &source = "") const = 0;

        // Add a new moving target to the database.
        //
        // The object must not already be defined.
        //
        // For const version, throws MovingTargetError if moving_target_id is
        // undefined.
        //
        // For non-const version, moving_target_id is updated if moving_target_id is
        // undefined.
        //
        // Throws `MovingTargetError` if this target or any of its names are
        // already in the database.
        virtual void add_moving_target(MovingTarget &target) const = 0;

        // Remove moving target from the database based on `moving_target_id`.
        virtual void remove_moving_target(const MovingTarget &target) const = 0;

        // Add a new observatory to the database that represents a particular
        // data source.
        //
        // Also updates the internal copy of the observatories, which is used
        // with `find_observations(Ephemeris)`;
        virtual void add_observatory(const string &name, const Observatory &observatory) const = 0;

        // Get an observatory from the database.
        virtual const Observatory get_observatory(const string &name) const = 0;

        // Get all observatories from the database.
        virtual const Observatories get_observatories() const = 0;

        // Remove an observatory from the database.
        //
        // Also updates the internal copy of the observatories, which is used
        // with `find_observations(Ephemeris)`;
        virtual void remove_observatory(const string &name) const = 0;

        // Get a list of observation sources.
        virtual const vector<string> get_sources() const = 0;

        // Update an existing moving target in the database based on `moving_target_id`.
        //
        // `moving_target_id` must be defined.
        //
        // Throws `MovingTargetNotFound` if the `moving_target_id` is not in the
        // database.
        virtual void update_moving_target(const MovingTarget &target) const = 0;

        // Get moving target by object ID.  Throws MovingTargetNotFound
        // if `moving_target_id` is not in database.
        virtual MovingTarget get_moving_target(const int moving_target_id) const = 0;

        // Get moving target by name and small body flag.  If name is not in the
        // database, returns a new MovingTarget object with an undefined
        // moving_target_id.
        virtual MovingTarget get_moving_target(const string &name, const bool small_body = true) const = 0;

        // Get all moving targets defined in the database.
        virtual vector<MovingTarget> get_all_moving_targets() const = 0;

        // Add ephemeris data to the database.
        //
        // The ephemeris's target ID must be in the database.
        virtual void add_ephemeris(Ephemeris &eph) const = 0;

        // Get ephemeris data from the database, optionally limited to a specific date range.
        virtual Ephemeris get_ephemeris(const MovingTarget target, double mjd_start = 0, double mjd_stop = 100000) const = 0;

        // Get the minimum and maximum dates of all ephemerides of all targets in the database.
        virtual std::pair<double *, double *> ephemeris_date_range() const = 0;

        // Remove ephemeris data from the database, optionally limited to a specific date range.
        virtual int remove_ephemeris(const MovingTarget target, double mjd_start = 0, double mjd_stop = 100000) const = 0;

        // Add an observation to the database.
        // - generally one would use sbsearch.add_observations()
        // - if observation ID is set, the database entry for this ID is updated
        // - if the observation ID is not set, a new database entry is made and
        //   the observation will be updated with the new ID
        // - index terms must be defined
        virtual void add_observation(Observation &observation) const = 0;

        // Add a set of observations to the database, see add_observation for details.
        void add_observations(Observations &observations) const;

        // Get an observation from the database.
        virtual Observation get_observation(const int64 observation_id) const = 0;

        // Get a set of observations from the database by observation_id, from first up to last.
        template <typename ForwardIterator>
        Observations get_observations(const ForwardIterator &first, const ForwardIterator &last) const;

        // Remove observations matching date range.
        virtual void remove_observations(const double mjd_start, const double mjd_stop) const = 0;

        // Remove observations matching source and date range.
        virtual void remove_observations(const string &source, const double mjd_start, const double mjd_stop) const = 0;

        // Search options.
        //
        // Found observations will be fully within the mjd limits.
        //
        // With parallax accounting enabled for ephemeris searches, the target
        // ephemeris must be computed for the geocenter, and the observatory
        // parallax constants defined.
        struct Options
        {
            double mjd_start = 0; // default: effectively search over all time
            double mjd_stop = 100000;
            string source = string(); // default: search all sources
            bool parallax = false;    // default: do not account for ephemeris parallax
        };

        // Count observations matching dates.
        virtual int64 count_observations(const double mjd_start, const double mjd_stop) const = 0;

        // Count observations matching source, and optionally dates.
        virtual int64 count_observations(const string &source, const double mjd_start, const double mjd_stop) const = 0;

        // Find observations by date.
        virtual Observations find_observations(const double mjd_start, const double mjd_stop, const int64 limit, const int64 offset) const = 0;

        // Find observations by source and date.
        virtual Observations find_observations(const string &source, const double mjd_start, double mjd_stop, const int64 limit, const int64 offset) const = 0;

        // Find observations matched by the provided query terms.
        virtual Observations find_observations(vector<string> query_terms, const Options &options) const = 0;

        // Add a found object to the database.
        virtual void add_found(const Found &found) const = 0;

        // Add found objects to the database.
        void add_founds(const Founds &founds) const;

        // Get all found moving targets for an observation from the database.
        virtual Founds get_found(const Observation &observation) const = 0;

        // Get all found observations for a moving target from the database.
        virtual Founds get_found(const MovingTarget &target) const = 0;

        // Remove a found object from the database.
        virtual void remove_found(const Found &found) const = 0;

        // Remove found objects from the database.
        void remove_founds(const Founds &founds) const;
    };

    template <typename ForwardIterator>
    Observations SBSearchDatabase::get_observations(const ForwardIterator &first_observation_id, const ForwardIterator &last_observation_id) const
    {
        Observations observations;
        ForwardIterator observation_id = first_observation_id;
        while (observation_id != last_observation_id)
        {
            observations.push_back(get_observation(*observation_id));
            observation_id++;
        }
        return observations;
    }
}
#endif // SBSDB_H_
