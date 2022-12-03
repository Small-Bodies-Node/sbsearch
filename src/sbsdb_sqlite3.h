#ifndef SBSDB_SQLITE3_H_
#define SBSDB_SQLITE3_H_

#include "ephemeris.h"
#include "observation.h"
#include "sbsdb.h"

#include "sqlite3.h"

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

        // initialize database, or add any missing tables, indices, etc.
        void setup_tables() override;

        // get a single value from a SQL statement
        template <typename T>
        T get_one_value(const char *statement);

        // add a single Observation to the database
        void add_observation(Observation obervation) override;

        // get a single Observation from the database
        Observation get_observation(const int64 observation_id) override;

        // using SBSearchDatabase::fuzzy_search;  // not working but why?
        vector<Observation> fuzzy_search(vector<string> terms) override;
        vector<Observation> fuzzy_search(Ephemeris eph) override;

    private:
        sqlite3 *db;
        void execute_sql(const char *statement) override;
        void check_sql(char *error_message);
    };

    // definte templates
    template <typename T>
    T SBSearchDatabaseSqlite3::get_one_value(const char *statement)
    {
        auto set_value = [](void *val, int count, char **data, char **columns)
        {
            T *converted_value = (T *)val; // convert void* to T*
            std::stringstream convert(data[0]);
            convert >> *converted_value;
            return 0;
        };

        char *error_message = NULL;
        T value;
        sqlite3_exec(db, statement, set_value, &value, &error_message);
        check_sql(error_message);
        return value;
    }

}
#endif // SBSDB_SQLITE3_H_