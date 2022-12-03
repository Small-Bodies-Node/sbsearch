#ifndef SBSDB_H_
#define SBSDB_H_

#include "ephemeris.h"
#include "observation.h"

#include <s2/s2point.h>
#include <s2/s2cap.h>
#include <s2/s2polyline.h>
#include <s2/s2polygon.h>
#include <s2/s2region_term_indexer.h>

namespace sbsearch
{
    class SBSearchDatabase
    {
    public:
        SBSearchDatabase();

        // initialize database, or add any missing tables, indices, etc.
        virtual void setup_tables() = 0;

        void drop_time_indices();

        // void add_moving_target(const MovingTarget target);
        // MovingTarget get_moving_target(const int64 object_id);
        // MovingTarget get_moving_target(const char* name);

        // add an ephemeris to the database
        // void add_ephemeris(const Ephemeris ephemeris);
        // get an ephemeris from the database
        // Ephemeris get_ephemeris(const MovingTarget target);

        // add an observation to the database, if the observation ID is not set, it will be updated
        void add_observation(Observation observation);
        // get an observation from the database
        Observation get_observation(const int64 observation_id);

        // add a set of observations to the database, they may be updated (see add_observation)
        void add_observations(vector<Observation> &observations);
        // get a set of observations from the database
        vector<Observation> get_observations(const vector<int64> observation_ids);

        vector<Observation> fuzzy_search(vector<string> terms);

        vector<Observation> find_observations(S2Point point);
        vector<Observation> find_observations(S2Cap cap);
        vector<Observation> find_observations(S2Polyline polyline);
        vector<Observation> find_observations(S2Polygon polygon);
        vector<Observation> find_observations(Ephemeris ephemeris);

    private:
        // this method does the actual database transaction
        virtual void add_observation_sql(Observation observation, const string terms_string) = 0;
        virtual void execute_sql(const char *statement) = 0;
        // void execute_sql(const char *statement);

        S2RegionTermIndexer indexer;
    };
}
#endif // SBSDB_H_
