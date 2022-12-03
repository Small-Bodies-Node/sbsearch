#include "sbsdb.h"
#include "sbsdb_sqlite3.h"
#include "observation.h"
#include "sbsearch.h"

#include <iostream>
#include <stdexcept>
#include <string>

#include "sqlite3.h"

using std::cerr;
using std::cout;
using std::endl;

// #define SQLERRCHECK                                    \
//     {                                                  \
//         if (error_message != NULL)                     \
//         {                                              \
//             printf("%s\n", error_message);             \
//             sqlite3_free(error_message);               \
//             throw std::runtime_error("\nSQL error\n"); \
//         }                                              \
//     }

namespace sbsearch
{
    SBSearchDatabaseSqlite3::SBSearchDatabaseSqlite3(const char *filename)
    {
        int rc = sqlite3_open(filename, &db);
        if (rc != SQLITE_OK)
        {
            cerr << "Error opening database" << sqlite3_errmsg(db) << endl;
            throw std::runtime_error("Error opening database");
        }
        else
            cout << "Opened " << filename << ".\n";
        execute_sql("PRAGMA temp_store_directory = './';");
    }

    SBSearchDatabaseSqlite3::~SBSearchDatabaseSqlite3()
    {
        sqlite3_close(db);
        cout << "Closed database.\n";
    }

    void SBSearchDatabaseSqlite3::setup_tables()
    {
        execute_sql("CREATE TABLE IF NOT EXISTS observations ("
                    "  observation_id INTEGER PRIMARY KEY,"
                    "  mjd_start FLOAT NOT NULL,"
                    "  mjd_stop FLOAT NOT NULL,"
                    "  fov TEXT NOT NULL,"
                    "  terms TEXT NOT NULL"
                    ");");
        execute_sql("CREATE VIRTUAL TABLE IF NOT EXISTS observations_geometry_time USING fts5(terms);");
        execute_sql("CREATE INDEX IF NOT EXISTS idx_observations_mjdstart ON observations(mjd_start);\n"
                    "CREATE INDEX IF NOT EXISTS idx_observations_mjdstop ON observations(mjd_stop);");
        execute_sql("CREATE TRIGGER IF NOT EXISTS on_insert_observations_insert_observations_geometry_time\n"
                    " AFTER INSERT ON observations\n"
                    " BEGIN\n"
                    "   INSERT INTO observations_geometry_time(rowid, terms)\n"
                    "   VALUES (new.observation_id, new.terms);\n"
                    " END;\n");
        execute_sql("CREATE TRIGGER IF NOT EXISTS on_update_observations_insert_observations_geometry_time\n"
                    " AFTER UPDATE OF observation_id, terms ON observations\n"
                    " BEGIN\n"
                    "   UPDATE observations_geometry_time\n"
                    "   SET rowid = new.observation_id, terms = new.terms\n"
                    "   WHERE rowid = old.observation_id;\n"
                    " END;");
        execute_sql("CREATE TRIGGER IF NOT EXISTS on_delete_observations_insert_observations_geometry_time\n"
                    " AFTER DELETE ON observations\n"
                    " BEGIN\n"
                    "   DELETE FROM observations_geometry_time\n"
                    "   WHERE observation_id = old.observation_id;\n"
                    " END;");

        cout << "Tables are set." << endl;
    }

    void SBSearchDatabaseSqlite3::check_sql(char *error_message)
    {
        if (error_message != NULL)
        {
            // printf("%s\n", error_message);
            cerr << error_message << endl;
            sqlite3_free(error_message);
            throw std::runtime_error("\nSQL error\n");
        }
    }

    void SBSearchDatabaseSqlite3::execute_sql(const char *statement)
    {
        char *error_message = NULL;
        sqlite3_exec(db, statement, NULL, 0, &error_message);
        check_sql(error_message);
    }

    void SBSearchDatabaseSqlite3::add_observation_sql(Observation observation, const string terms)
    {
        char *error_message = NULL;
        int rc;
        int64 observation_id;
        sqlite3_stmt *statement;

        if (observation.observation_id() == -1)
        {
            // insert row and update observation object with observation_id
            rc = sqlite3_prepare_v2(db, "INSERT INTO observations VALUES (NULL, ?, ?, ?, ?) RETURNING observation_id;", -1, &statement, NULL);
        }
        else
        {
            // update existing observation
            rc = sqlite3_prepare_v2(db, "UPDATE observations SET mjd_start=?, mjd_stop=?, fov=?, terms=? WHERE observation_id=?;", -1, &statement, NULL);
            rc += sqlite3_bind_int64(statement, 5, observation.observation_id());
        }

        rc += sqlite3_bind_double(statement, 1, observation.mjd_start());
        rc += sqlite3_bind_double(statement, 2, observation.mjd_stop());
        rc += sqlite3_bind_text(statement, 3, observation.fov().c_str(), -1, SQLITE_TRANSIENT);
        rc += sqlite3_bind_text(statement, 4, terms.c_str(), -1, SQLITE_TRANSIENT);
        check_sql(error_message);

        rc = sqlite3_step(statement);
        check_sql(error_message);

        if (observation.observation_id() == -1)
            observation.observation_id(sqlite3_column_int64(statement, 0));
        else if (rc != SQLITE_DONE)
        {
            std::cerr << sqlite3_errmsg(db) << std::endl;
            throw std::runtime_error("Error updating observation in database");
        }

        sqlite3_finalize(statement);

        // rc = sqlite3_prepare_v2(db, "INSERT INTO obs_geometry_time(rowid, terms) VALUES (?, ?);", -1, &statement, NULL);
        // rc += sqlite3_bind_int64(statement, 1, observation.observation_id());
        // rc += sqlite3_bind_text(statement, 3, terms.c_str(), -1, SQLITE_TRANSIENT);
        // if (rc != SQLITE_OK)
        // {
        //     cerr << sqlite3_errmsg(db) << endl;
        //     throw std::runtime_error("Error preparing SQL statement");
        // }

        // rc = sqlite3_step(statement);
        // if (rc != SQLITE_DONE)
        // {
        //     std::cerr << sqlite3_errmsg(db) << std::endl;
        //     throw std::runtime_error("Error updating observation in database");
        // }

        // sqlite3_finalize(statement);
    }
}