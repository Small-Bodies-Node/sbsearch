#ifndef SBSEARCH_H_
#define SBSEARCH_H_

#include <string>
#include <vector>
#include <s2/s2latlng.h>
#include <s2/s2region.h>

#include "ephemeris.h"
#include "indexer.h"
#include "sbsdb.h"

namespace sbsearch
{
    struct Found
    {
        Observation observation;
        Ephemeris ephemeris;

        Found(Observation o, Ephemeris e) : observation(o), ephemeris(e){};
    };

    std::ostream &operator<<(std::ostream &os, const Found &found);
    // If show_fov is set on any observation, then the FOV is shown for all
    std::ostream &operator<<(std::ostream &os, const vector<Found> &founds);

    class SBSearch
    {
    public:
        enum DatabaseType
        {
            sqlite3
        };

        // constructor
        //
        // Setting log_file has no effect if the Logger has already been initalized.
        //
        // For sqlite3 databases:
        //   - `name` is the database filename, ":memory:" for an in-memory
        //     database, or "" (empty-string) for a temporary on-disk database.
        SBSearch(DatabaseType database_type, const std::string name, const std::string log_file = "/dev/null");

        // database maintainence
        //
        // drop/create indices, generally used when adding many new observations
        inline void drop_observations_indices() { db->drop_observations_indices(); };
        inline void create_observations_indices() { db->create_observations_indices(); };

        // read-only access to indexer options
        const Indexer::Options &indexer_options() { return indexer_.options(); };

        // Re-index the terms for each observation and ephemeris, and
        // store the new indexer parameters to the database.
        void reindex(Indexer::Options options);

        // database I/O

        // Add moving target, object_id will be updated as needed.
        void add_moving_target(MovingTarget &target);

        // Remove moving target from the database based on `object_id`.
        void remove_moving_target(const MovingTarget &target);

        // Update an existing moving target in the database based on `object_id`.
        //
        // `object_id` must be defined.
        void update_moving_target(const MovingTarget &target);

        // Get moving target by object ID or name.
        MovingTarget get_moving_target(const int object_id);
        MovingTarget get_moving_target(const string &name);

        // Add ephemeris data to the database.
        //
        // If the ephemeris's target is not already in the database, then it
        // will be added and eph.target() updated.
        void add_ephemeris(Ephemeris &eph);

        // Get ephemeris data from the database, optionally limited to a specific date range.
        Ephemeris get_ephemeris(const MovingTarget target, double mjd_start = 0, double mjd_stop = 70000);

        // Remove ephemeris data from the database, optionally limited to a specific date range.
        int remove_ephemeris(const MovingTarget target, double mjd_start = 0, double mjd_stop = 70000);

        // Add observations, index terms will be added as needed.
        void add_observations(vector<Observation> &observations);
        vector<Observation> get_observations(const vector<int64> &observation_id);

        // Start and end dates, optionally for a specific survey.
        std::pair<double *, double *> date_range(string source = "");

        // search functions

        // Search by point or polygon, optionally over a time range.  Let
        // `mjd_start` or `mjd_stop` = -1 for an unbounded limit.
        vector<Observation> find_observations(const S2Point &point, double mjd_start = -1, double mjd_stop = -1);
        vector<Observation> find_observations(const S2Polygon &polygon, double mjd_start = -1, double mjd_stop = -1);

        // Search by ephemeris.
        vector<Found> find_observations(const Ephemeris &ephemeris);

    private:
        SBSearchDatabase *db;
        Indexer indexer_;
    };
}

#endif // SBSEARCH_H_
