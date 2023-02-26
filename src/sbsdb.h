#ifndef SBSDB_H_
#define SBSDB_H_

#include <s2/s2point.h>
#include <s2/s2cap.h>
#include <s2/s2metrics.h>
#include <s2/s2polyline.h>
#include <s2/s2polygon.h>
#include <s2/s2region_term_indexer.h>

#include "ephemeris.h"
#include "indexer.h"
#include "observation.h"

#define SBSEARCH_DATABASE_VERSION "3.0"

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
        virtual double get_double(const char *statement) = 0;
        virtual int get_int(const char *statement) = 0;
        virtual int64 get_int64(const char *statement) = 0;
        virtual string get_string(const char *statement) = 0;

        // drop/create observations indices, e.g., before inserting many new observations
        // restore indices with setup_tables()
        void drop_observations_indices();
        void create_observations_indices();

        // Read indexer options from the database.
        Indexer::Options indexer_options();
        // Write indexer options to the database.
        virtual void indexer_options(Indexer::Options options) = 0;

        // get date range, optionally for a single source
        virtual std::pair<double, double> date_range(string source = "") = 0;

        // void add_moving_target(const MovingTarget target);
        // MovingTarget get_moving_target(const int64 object_id);
        // MovingTarget get_moving_target(const char* name);

        // add an ephemeris to the database
        // void add_ephemeris(const Ephemeris ephemeris);
        // get an ephemeris from the database
        // Ephemeris get_ephemeris(const MovingTarget target);

        // add an observation to the database
        // - if the observation ID is not set, it will be updated
        // - index terms must be defined
        virtual void add_observation(Observation &observation) = 0;

        // add a set of observations to the database, see add_observation for details
        void add_observations(vector<Observation> &observations);

        // get an observation from the database
        virtual Observation get_observation(const int64 observation_id) = 0;

        // get a set of observations from the database by observation_id, from first up to last.
        template <typename ForwardIterator>
        vector<Observation> get_observations(const ForwardIterator &first, const ForwardIterator &last);

        // Find observations matched by the provided query terms.
        virtual vector<Observation> find_observations(vector<string> query_terms) = 0;
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
