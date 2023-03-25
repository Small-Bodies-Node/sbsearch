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
// 10 to 10000... 10000 was slowest, but the difference is really small
#define MAXIMUM_QUERY_TERMS size_t(100)

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

        void add_observatory(const string &name, const Observatory &observatory) override;
        const Observatory get_observatory(const string &name) override;
        const Observatories get_observatories() override;
        void remove_observatory(const string &name) override;

        void add_ephemeris(Ephemeris &eph) override;
        Ephemeris get_ephemeris(const MovingTarget target, double mjd_start = 0, double mjd_stop = 70000) override;
        int remove_ephemeris(const MovingTarget target, double mjd_start = 0, double mjd_stop = 70000) override;

        void add_observation(Observation &observation) override;
        Observation get_observation(const int64 observation_id) override;
        vector<Observation> find_observations(vector<string> query_terms, const Options &options = Options()) override;

    private:
        sqlite3 *db;
        void check_rc(const int rc);
        void check_sql(char *error_message);
        void error_if_closed();
        void add_moving_target_name(const int object_id, const string &name, const bool primary_id);
    };
}
#endif // SBSDB_SQLITE3_H_