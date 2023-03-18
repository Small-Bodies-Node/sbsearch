#ifndef SBSDB_SQLITE3_H_
#define SBSDB_SQLITE3_H_

#include "observation.h"
#include "sbsdb.h"

#include <string>
#include <sqlite3.h>
#include <s2/s2point.h>
#include <s2/s2cap.h>
#include <s2/s2polyline.h>
#include <s2/s2polygon.h>
#include <s2/s2region_term_indexer.h>

// Order 1000 seems fine, 10,000 was slower in testing
#define MAXIMUM_QUERY_CLAUSE_LENGTH 1000

namespace sbsearch
{
    class SBSearchDatabaseSqlite3 : public SBSearchDatabase
    {
    public:
        SBSearchDatabaseSqlite3(const std::string filename);
        ~SBSearchDatabaseSqlite3()
        {
            close();
        }

        void close() override;

        void setup_tables() override;

        void execute_sql(const char *statement) override;

        // sqlite's execute with callback
        void execute_sql(const char *statement, int (*callback)(void *, int, char **, char **), void *callback_arg);

        double *get_double(const char *statement) override;
        int *get_int(const char *statement) override;
        int64 *get_int64(const char *statement) override;
        string *get_string(const char *statement) override;

        void indexer_options(Indexer::Options options) override;

        std::pair<double *, double *> date_range(const string &source = "") override;

        void add_moving_target(MovingTarget &target) override;
        void remove_moving_target(const MovingTarget &target) override;
        void update_moving_target(const MovingTarget &target) override;
        MovingTarget get_moving_target(const int object_id) override;
        MovingTarget get_moving_target(const string &name) override;

        // add an observation to the database
        // - generally one would use sbsearch.add_observations()
        // - if observation ID is set, the database entry for this ID is updated
        // - if the observation ID is not set, a new database entry is made and
        //   the observation will be updated with the new ID
        // - index terms must be defined
        void add_observation(Observation &observation) override;

        // get a single Observation from the database
        Observation get_observation(const int64 observation_id) override;

        vector<Observation> find_observations(vector<string> query_terms) override;

    private:
        sqlite3 *db;
        void check_rc(const int rc);
        void check_sql(char *error_message);
        void error_if_closed();
        void add_moving_target_name(const int object_id, const string &name, const bool primary_id);
    };
}
#endif // SBSDB_SQLITE3_H_