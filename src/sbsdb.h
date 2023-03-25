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
        virtual void execute_sql(const char *statement) = 0;

        // get single value results from a SQL statement
        virtual double *get_double(const char *statement) = 0;
        virtual int *get_int(const char *statement) = 0;
        virtual int64 *get_int64(const char *statement) = 0;
        virtual string *get_string(const char *statement) = 0;

        // drop/create observations indices, e.g., before inserting many new observations
        // restore indices with setup_tables()
        void drop_observations_indices();
        void create_observations_indices();

        // Read indexer options from the database.
        Indexer::Options indexer_options();
        // Write indexer options to the database.
        virtual void indexer_options(Indexer::Options options) = 0;

        // get date range, optionally for a single source
        virtual std::pair<double *, double *> date_range(const string &source = "") = 0;

        // Add a new moving target to the database.
        //
        // The object must not already be defined.
        //
        // For const version, throws MovingTargetError if object_id is
        // undefined.
        //
        // For non-const version, object_id is updated if object_id is
        // undefined.
        //
        // Throws `MovingTargetError` if this target or any of its names are
        // already in the database.
        virtual void add_moving_target(MovingTarget &target) = 0;

        // Remove moving target from the database based on `object_id`.
        virtual void remove_moving_target(const MovingTarget &target) = 0;

        // Add a new observatory to the database that represents a particular data source.
        virtual void add_observatory(const string &source, const Observatory &observatory) = 0;

        // Get an observatory from the database.
        virtual const Observatory get_observatory(const string &source) = 0;

        // Get all observatories from the database.
        virtual const Observatories get_observatories() = 0;

        // Remove an observatory from the database.
        virtual void remove_observatory(const string &source) = 0;

        // Update an existing moving target in the database based on `object_id`.
        //
        // `object_id` must be defined.
        //
        // Throws `MovingTargetNotFound` if the `object_id` is not in the
        // database.
        virtual void update_moving_target(const MovingTarget &target) = 0;

        // Get moving target by object ID or name.
        //
        // Throws MovingTargetNotFound if `object_id` or `name` is not in database.
        virtual MovingTarget get_moving_target(const int object_id) = 0;
        virtual MovingTarget get_moving_target(const string &name) = 0;

        // Add ephemeris data to the database.
        //
        // If the ephemeris's target is not already in the database, then it
        // will be added and eph.target() updated.
        virtual void add_ephemeris(Ephemeris &eph) = 0;

        // Get ephemeris data from the database, optionally limited to a specific date range.
        virtual Ephemeris get_ephemeris(const MovingTarget target, double mjd_start = 0, double mjd_stop = 70000) = 0;

        // Remove ephemeris data from the database, optionally limited to a specific date range.
        virtual int remove_ephemeris(const MovingTarget target, double mjd_start = 0, double mjd_stop = 70000) = 0;

        // Add an observation to the database.
        // - generally one would use sbsearch.add_observations()
        // - if observation ID is set, the database entry for this ID is updated
        // - if the observation ID is not set, a new database entry is made and
        //   the observation will be updated with the new ID
        // - index terms must be defined
        virtual void add_observation(Observation &observation) = 0;

        // Add a set of observations to the database, see add_observation for details.
        void add_observations(vector<Observation> &observations);

        // Get an observation from the database.
        virtual Observation get_observation(const int64 observation_id) = 0;

        // Get a set of observations from the database by observation_id, from first up to last.
        template <typename ForwardIterator>
        vector<Observation> get_observations(const ForwardIterator &first, const ForwardIterator &last);

        // Search options.
        //
        // Observations must be fully within the mjd limits.
        //
        // With parallax accounting enabled for ephemeris searches, the target
        // must be computed for the geocenter, and the observatory parallax
        // constants defined.
        struct Options
        {
            double mjd_start = 0; // default: effectively search over all time
            double mjd_stop = 100000;
            string source = string();    // default: search all sources
            bool parallax = false;       // default: do not account for ephemeris parallax
            Observatories observatories; // parallax requires observatories keyed by source name
        };

        // Find observations matched by the provided query terms.
        virtual vector<Observation>
        find_observations(vector<string> query_terms, const Options &options) = 0;
    };

    template <typename ForwardIterator>
    vector<Observation> SBSearchDatabase::get_observations(const ForwardIterator &first_observation_id, const ForwardIterator &last_observation_id)
    {
        vector<Observation> observations;
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
