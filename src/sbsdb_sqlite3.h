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

namespace sbsearch
{
    class SBSearchDatabaseSqlite3 : public SBSearchDatabase
    {
    public:
        SBSearchDatabaseSqlite3(const char *filename);
        ~SBSearchDatabaseSqlite3();

        // initialize database, or add any missing tables, indices, etc.
        void setup_tables() override;

    private:
        sqlite3 *db;
        void add_observation_sql(Observation observation, const string terms) override;
        void execute_sql(const char *statement) override;
        void check_sql(char *error_message);
    };
}

#endif // SBSDB_SQLITE3_H_