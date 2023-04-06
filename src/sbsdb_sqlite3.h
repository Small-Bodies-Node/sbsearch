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

        void execute_sql(const char *statement) const override;

        // sqlite's execute with callback
        void execute_sql(const char *statement, int (*callback)(void *, int, char **, char **), void *callback_arg) const;

        double *get_double(const char *statement) const override;
        int *get_int(const char *statement) const override;
        int64 *get_int64(const char *statement) const override;
        string *get_string(const char *statement) const override;

        void indexer_options(Indexer::Options options) override;

        std::pair<double *, double *> observation_date_range(const string &source = "") const override;

        void add_moving_target(MovingTarget &target) const override;
        void remove_moving_target(const MovingTarget &target) const override;
        void update_moving_target(const MovingTarget &target) const override;
        MovingTarget get_moving_target(const int moving_target_id) const override;
        MovingTarget get_moving_target(const string &name) const override;
        vector<MovingTarget> get_all_moving_targets() const override;

        void add_observatory(const string &name, const Observatory &observatory) const override;
        const Observatory get_observatory(const string &name) const override;
        const Observatories get_observatories() const override;
        void remove_observatory(const string &name) const override;
        const vector<string> get_sources() const override;

        void add_ephemeris(Ephemeris &eph) const override;
        Ephemeris get_ephemeris(const MovingTarget target, double mjd_start = 0, double mjd_stop = 100000) const override;
        int remove_ephemeris(const MovingTarget target, double mjd_start = 0, double mjd_stop = 100000) const override;
        std::pair<double *, double *> ephemeris_date_range() const override;

        void add_observation(Observation &observation) const override;
        Observation get_observation(const int64 observation_id) const override;
        void remove_observations(const double mjd_start, const double mjd_stop) const override;
        void remove_observations(const string &source, const double mjd_start, const double mjd_stop) const override;

        int64 count_observations(const double mjd_start, const double mjd_stop) const override;
        int64 count_observations(const string &source, const double mjd_start, const double mjd_stop) const override;

        Observations find_observations(const double mjd_start, const double mjd_stop, const int64 limit, const int64 offset) const override;
        Observations find_observations(const string &source, const double mjd_start, double mjd_stop, const int64 limit, const int64 offset) const override;
        Observations find_observations(vector<string> query_terms, const Options &options = Options()) const override;

        void add_found(const Found &found) const override;
        vector<Found> get_found(const Observation &observation) const override;
        vector<Found> get_found(const MovingTarget &target) const override;
        void remove_found(const Found &found) const override;

    private:
        sqlite3 *db;
        void check_rc(const int rc) const;
        void check_sql(char *error_message) const;
        void error_if_closed() const;
        void add_moving_target_name(const int moving_target_id, const string &name, const bool primary_id) const;
    };
}
#endif // SBSDB_SQLITE3_H_