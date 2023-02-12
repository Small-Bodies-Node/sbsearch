#ifndef SBSDB_SQLITE3_H_
#define SBSDB_SQLITE3_H_

#include "ephemeris.h"
#include "observation.h"
#include "sbsdb.h"

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
        SBSearchDatabaseSqlite3(const char *filename);
        ~SBSearchDatabaseSqlite3();

        // close the database connection
        void close() override;

        // initialize database, or add any missing tables, indices, etc.
        void setup_tables() override;

        // execute a sql statement
        void execute_sql(const char *statement) override;
        void execute_sql(const char *statement, int (*callback)(void *, int, char **, char **), void *callback_arg);

        // get a single string result from a SQL statement
        double get_double(const char *statement) override;
        int get_int(const char *statement) override;
        int64 get_int64(const char *statement) override;
        string get_string(const char *statement) override;

        // get date range, optionally for a single source
        std::pair<double, double> date_range(string source = "") override;

        // add an observation to the database
        // - generally one would use sbsearch.add_observations()
        // - if the observation ID is not set, it will be updated
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
        template <typename T>
        static int get_single_value_callback(void *val, int count, char **data, char **columns);
    };

    // templated functions
    template <typename T>
    int SBSearchDatabaseSqlite3::get_single_value_callback(void *val, int count, char **data, char **columns)
    {
        T *converted_value = (T *)val;
        std::stringstream convert(data[0]);
        convert >> *converted_value;
        return 0;
    };
}
#endif // SBSDB_SQLITE3_H_