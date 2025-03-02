#ifndef SBSDB_POSTGRESQL_H_
#define SBSDB_POSTGRESQL_H_

#include "observation.h"
#include "sbsdb.h"

#include <string>
#include <pqxx/pqxx>
#include <sqlite3.h>
#include <s2/s2point.h>
#include <s2/s2cap.h>
#include <s2/s2polyline.h>
#include <s2/s2polygon.h>
#include <s2/s2region_term_indexer.h>

// Order 1000?  Need to test
#define MAXIMUM_QUERY_CLAUSE_LENGTH 1000
// 10 to 10000... Need to test
#define MAXIMUM_QUERY_TERMS size_t(100)

namespace sbsearch
{
    class SBSearchDatabasePostgreSQL : public SBSearchDatabase
    {
    public:
        SBSearchDatabasePostgreSQL(const std::string &uri) : connection_(uri)
        {
            if (connection_.is_open())
                Logger::info() << "Opened database URI " << uri << std::endl;
            else
                throw "Failed to open database connection.";
        }

        ~SBSearchDatabasePostgreSQL()
        {
            close();
        }

        void close() override;

        void setup_tables() override;

        void execute_sql(const char *statement) override;

        void drop_observations_indices() override;
        void create_observations_indices() override;

        // get single value results from a SQL statement
        double *get_double(const char *statement) override;
        int *get_int(const char *statement) override;
        int64 *get_int64(const char *statement) override;
        string *get_string(const char *statement) override;

        void indexer_options(Indexer::Options options) override;

        std::pair<double *, double *> observation_date_range(const string &source = "") override;

        void add_moving_target(MovingTarget &target) override;
        void remove_moving_target(const MovingTarget &target) override;
        MovingTarget get_moving_target(const int moving_target_id) override;
        MovingTarget get_moving_target(const string &name, const bool small_body = true) override;
        vector<MovingTarget> get_all_moving_targets() override;

        void add_observatory(const string &name, const Observatory &observatory) override;
        const Observatory get_observatory(const string &name) override;
        const Observatories get_observatories() override;
        void remove_observatory(const string &name) override;
        const vector<string> get_sources() override;

        void add_ephemeris(Ephemeris &eph) override;
        Ephemeris get_ephemeris(const MovingTarget target, double mjd_start = 0, double mjd_stop = 100000) override;
        int remove_ephemeris(const MovingTarget target, double mjd_start = 0, double mjd_stop = 100000) override;

        void add_observation(Observation &observation) override;
        Observation get_observation(const int64 observation_id) override;
        void remove_observations(const double mjd_start, const double mjd_stop) override;
        void remove_observations(const string &source, const double mjd_start, const double mjd_stop) override;

        // Count number of observations within an interval.
        int64 count_observations(const double mjd_start, const double mjd_stop) override;

        // Count number of observations for a source within an interval, if
        // source is an empty string, then count all sources.
        int64 count_observations(const string &source, const double mjd_start, const double mjd_stop) override;

        Observations find_observations(const double mjd_start, const double mjd_stop, const int64 limit, const int64 offset) override;
        Observations find_observations(const string &source, const double mjd_start, double mjd_stop, const int64 limit, const int64 offset) override;
        Observations find_observations(vector<string> query_terms, const Options &options = Options()) override;

        void add_found(const Found &found) override;
        Founds get_found(const Observation &observation) override;
        Founds get_found(const MovingTarget &target) override;
        void remove_found(const Found &found) override;

    private:
        pqxx::connection connection_;
        // pqxx::connection *connection_;
        void error_if_closed();
        void add_moving_target_name(pqxx::transaction_base &work, const int moving_target_id, const string &name, const bool small_body, const bool primary_id);
    };
}
#endif // SBSDB_POSTGRESQL_H_