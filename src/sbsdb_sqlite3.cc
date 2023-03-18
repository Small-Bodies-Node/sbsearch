#include "config.h"

#include <iostream>
#include <stdexcept>
#include <set>
#include <string>
#include <vector>
#include <sqlite3.h>

#include "logging.h"
#include "observation.h"
#include "sbsdb.h"
#include "sbsdb_sqlite3.h"

using std::endl;
using std::set;
using std::string;
using std::vector;

namespace sbsearch
{
    SBSearchDatabaseSqlite3::SBSearchDatabaseSqlite3(const string filename)
    {
        int rc = sqlite3_open(filename.c_str(), &db);
        if (rc != SQLITE_OK)
        {
            Logger::error() << "Error opening database: " << sqlite3_errmsg(db) << endl;
            throw std::runtime_error("Error opening database");
        }
        else
            Logger::info() << "Opened sqlite3 database " << filename << endl;
        execute_sql("PRAGMA temp_store_directory = './';");
    }

    void SBSearchDatabaseSqlite3::close()
    {
        if (db != NULL)
        {
            sqlite3_close(db);
            db = NULL;
            Logger::info() << "Closed database." << endl;
        }
    }

    void SBSearchDatabaseSqlite3::setup_tables()
    {
        execute_sql(R"(
CREATE TABLE IF NOT EXISTS observations (
  observation_id INTEGER PRIMARY KEY,
  source TEXT NOT NULL,
  product_id TEXT NOT NULL,
  mjd_start FLOAT NOT NULL,
  mjd_stop FLOAT NOT NULL,
  fov TEXT NOT NULL,
  terms TEXT NOT NULL
);
)");
        execute_sql("CREATE VIRTUAL TABLE IF NOT EXISTS observations_terms_index USING fts5(terms, content='observations', content_rowid='observation_id');");

        // Triggers to keep the FTS index up to date.
        execute_sql(R"(
CREATE TRIGGER IF NOT EXISTS observations_insert AFTER INSERT ON observations BEGIN
  INSERT INTO observations_terms_index(rowid, terms) VALUES (new.observation_id, new.terms);
END;
CREATE TRIGGER IF NOT EXISTS observations_delete AFTER DELETE ON observations BEGIN
  INSERT INTO observations_terms_index(observations_terms_index, rowid, terms) VALUES('delete', old.observation_id, old.terms);
END;
CREATE TRIGGER IF NOT EXISTS observations_update AFTER UPDATE ON observations BEGIN
  INSERT INTO observations_terms_index(observations_terms_index, rowid, terms) VALUES('delete', old.observation_id, old.terms);
  INSERT INTO observations_terms_index(rowid, terms) VALUES (new.observation_id, new.terms);
END;
)");

        execute_sql(R"(
CREATE TABLE IF NOT EXISTS moving_targets (
  moving_target_id INTEGER PRIMARY KEY,
  object_id INTEGER NOT NULL,
  name TEXT UNIQUE NOT NULL,
  primary_id BOOLEAN NOT NULL
);
CREATE UNIQUE INDEX IF NOT EXISTS moving_target_primary_id ON moving_targets(object_id, name) WHERE primary_id=TRUE;
CREATE INDEX IF NOT EXISTS moving_target_object_id ON moving_targets(object_id);
)");

        create_observations_indices();

        // configuration table and defaults
        execute_sql(R"(
CREATE TABLE IF NOT EXISTS configuration (
  parameter TEXT NOT NULL UNIQUE,
  value TEXT NOT NULL
);
INSERT OR IGNORE INTO configuration VALUES ('max_spatial_cells', '8');
INSERT OR IGNORE INTO configuration VALUES ('max_spatial_level', '12');
INSERT OR IGNORE INTO configuration VALUES ('min_spatial_level', '4');
INSERT OR IGNORE INTO configuration VALUES ('temporal_resolution', '1');
INSERT OR IGNORE INTO configuration VALUES ('database version', ')" SBSEARCH_DATABASE_VERSION "');");

        Logger::debug() << "Database tables are set." << endl;
    }

    void SBSearchDatabaseSqlite3::execute_sql(const char *statement)
    {
        execute_sql(statement, NULL, NULL);
    }

    void SBSearchDatabaseSqlite3::execute_sql(const char *statement, int (*callback)(void *, int, char **, char **), void *callback_arg)
    {
        error_if_closed();

        char *error_message = NULL;
        int rc = sqlite3_exec(db, statement, callback, callback_arg, &error_message);
        check_sql(error_message);
    }

    double *SBSearchDatabaseSqlite3::get_double(const char *statement)
    {
        double *value = new double;

        sqlite3_stmt *stmt;
        int rc;
        rc = sqlite3_prepare_v2(db, statement, -1, &stmt, NULL);
        check_rc(rc);

        rc = sqlite3_step(stmt);
        check_rc(rc);

        if (sqlite3_column_type(stmt, 0) == SQLITE_NULL)
            value = nullptr;
        else
            *value = sqlite3_column_double(stmt, 0);

        sqlite3_finalize(stmt);

        return std::move(value);
    }

    int *SBSearchDatabaseSqlite3::get_int(const char *statement)
    {
        int *value = new int;

        sqlite3_stmt *stmt;
        int rc;
        rc = sqlite3_prepare_v2(db, statement, -1, &stmt, NULL);
        check_rc(rc);

        rc = sqlite3_step(stmt);
        check_rc(rc);

        if (sqlite3_column_type(stmt, 0) == SQLITE_NULL)
            value = nullptr;
        else
            *value = sqlite3_column_int(stmt, 0);

        sqlite3_finalize(stmt);

        return value;
    }

    int64 *SBSearchDatabaseSqlite3::get_int64(const char *statement)
    {
        int64 *value = new int64;

        sqlite3_stmt *stmt;
        int rc;
        rc = sqlite3_prepare_v2(db, statement, -1, &stmt, NULL);
        check_rc(rc);

        rc = sqlite3_step(stmt);
        check_rc(rc);

        if (sqlite3_column_type(stmt, 0) == SQLITE_NULL)
            value = nullptr;
        else
            *value = sqlite3_column_int64(stmt, 0);

        sqlite3_finalize(stmt);

        return value;
    }

    string *SBSearchDatabaseSqlite3::get_string(const char *statement)
    {
        string *value = new string();

        sqlite3_stmt *stmt;
        int rc;
        rc = sqlite3_prepare_v2(db, statement, -1, &stmt, NULL);
        check_rc(rc);

        rc = sqlite3_step(stmt);
        check_rc(rc);

        if (sqlite3_column_type(stmt, 0) == SQLITE_NULL)
            value = nullptr;
        else
            *value = string(reinterpret_cast<const char *>(sqlite3_column_text(stmt, 0)));

        sqlite3_finalize(stmt);

        return value;
    }

    void SBSearchDatabaseSqlite3::indexer_options(Indexer::Options options)
    {
        error_if_closed();

        int rc;
        sqlite3_stmt *statement;

        vector<string> parameters = {"max_spatial_cells",
                                     "max_spatial_level",
                                     "min_spatial_level",
                                     "temporal_resolutionstd::"};
        vector<string> values = {std::to_string(options.max_spatial_cells()),
                                 std::to_string(options.max_spatial_level()),
                                 std::to_string(options.min_spatial_level()),
                                 std::to_string(options.temporal_resolution())};
        for (int i = 0; i < parameters.size(); i++)
        {
            sqlite3_prepare_v2(db, "UPDATE configuration SET parameter=?1, value=?2 WHERE parameter=?1;", -1, &statement, NULL);
            sqlite3_bind_text(statement, 1, parameters[i].c_str(), parameters[i].size(), SQLITE_STATIC);
            sqlite3_bind_text(statement, 2, values[i].c_str(), values[i].size(), SQLITE_STATIC);
            rc = sqlite3_step(statement);
            sqlite3_finalize(statement);
        }
    }

    std::pair<double *, double *> SBSearchDatabaseSqlite3::date_range(const string &source)
    {
        double *mjd_start = new double;
        double *mjd_stop = new double;
        if (source == "")
        {
            mjd_start = get_double("SELECT MIN(mjd_start) FROM observations;");
            mjd_stop = get_double("SELECT MAX(mjd_stop) FROM observations;");
        }
        else
        {
            sqlite3_stmt *statement;
            int rc;
            rc = sqlite3_prepare_v2(db, "SELECT MIN(mjd_start), MAX(mjd_stop) FROM observations WHERE source=?;", -1, &statement, NULL);
            check_rc(rc);
            rc = sqlite3_bind_text(statement, 1, source.c_str(), -1, SQLITE_TRANSIENT);
            check_rc(rc);
            rc = sqlite3_step(statement);
            check_rc(rc);

            if (sqlite3_column_type(statement, 0) == SQLITE_NULL)
                mjd_start = nullptr;
            else
                *mjd_start = sqlite3_column_double(statement, 0);

            if (sqlite3_column_type(statement, 1) == SQLITE_NULL)
                mjd_stop = nullptr;
            else
                *mjd_stop = sqlite3_column_double(statement, 1);

            sqlite3_finalize(statement);
        }

        return std::pair<double *, double *>(std::move(mjd_start), std::move(mjd_stop));
    }

    void SBSearchDatabaseSqlite3::add_moving_target(MovingTarget &target)
    {
        error_if_closed();

        int object_id = target.object_id();
        if (object_id == UNDEF_OBJECT_ID)
        {
            // new objects use db largest object_id + 1, or 1 if there are no objects
            object_id = *get_int("SELECT IFNULL(MAX(object_id), 0) + 1 FROM moving_targets");
            target.object_id(object_id);
        }

        int rc;
        sqlite3_stmt *stmt;
        sqlite3_prepare_v2(db, "SELECT COUNT() FROM moving_targets WHERE object_id=?", -1, &stmt, NULL);
        sqlite3_bind_int(stmt, 1, target.object_id());
        rc = sqlite3_step(stmt);
        check_rc(rc);
        int count = sqlite3_column_int(stmt, 0);
        if (count != 0)
            throw MovingTargetError("object_id " + std::to_string(target.object_id()) + " already exists");
        sqlite3_finalize(stmt);

        Logger::info() << "Add moving target " << target.designation() << std::endl;
        execute_sql("BEGIN TRANSACTION;");
        add_moving_target_name(object_id, target.designation(), true);
        for (const string &name : target.alternate_names())
            add_moving_target_name(target.object_id(), name, false);
        execute_sql("END TRANSACTION;");
    }

    void SBSearchDatabaseSqlite3::add_moving_target_name(const int object_id, const string &name, const bool primary_id)
    {
        error_if_closed();
        int rc;
        sqlite3_stmt *stmt;
        sqlite3_prepare_v2(db, "INSERT INTO moving_targets (object_id, name, primary_id) VALUES (?, ?, ?)", -1, &stmt, NULL);
        sqlite3_bind_int(stmt, 1, object_id);
        sqlite3_bind_text(stmt, 2, name.c_str(), -1, SQLITE_TRANSIENT);
        sqlite3_bind_int(stmt, 3, primary_id);
        rc = sqlite3_step(stmt);
        check_rc(rc);
        sqlite3_finalize(stmt);
    };

    void SBSearchDatabaseSqlite3::remove_moving_target(const MovingTarget &target)
    {
        error_if_closed();
        int rc;
        sqlite3_stmt *stmt;

        execute_sql("BEGIN TRANSACTION;");
        sqlite3_prepare_v2(db, "SELECT COUNT() FROM moving_targets WHERE object_id=?", -1, &stmt, NULL);
        sqlite3_bind_int(stmt, 1, target.object_id());
        rc = sqlite3_step(stmt);
        check_rc(rc);
        int count = sqlite3_column_int(stmt, 0);
        if (count == 0)
            throw MovingTargetError("object_id " + std::to_string(target.object_id()) + " not found");
        sqlite3_finalize(stmt);

        Logger::info() << "Remove moving target " << target.designation() << std::endl;
        sqlite3_prepare_v2(db, "DELETE FROM moving_targets WHERE object_id=?", -1, &stmt, NULL);
        sqlite3_bind_int(stmt, 1, target.object_id());
        rc = sqlite3_step(stmt);
        check_rc(rc);
        sqlite3_finalize(stmt);

        execute_sql("END TRANSACTION;");
    }

    void SBSearchDatabaseSqlite3::update_moving_target(const MovingTarget &target)
    {
        Logger::info() << "Update moving target " << target.designation() << std::endl;
        remove_moving_target(target);
        MovingTarget copy(target);
        add_moving_target(copy);
    }

    MovingTarget SBSearchDatabaseSqlite3::get_moving_target(const int object_id)
    {
        error_if_closed();
        MovingTarget target;

        sqlite3_stmt *stmt;
        int rc;

        sqlite3_prepare_v2(db, "SELECT name, primary_id, object_id FROM moving_targets WHERE object_id=?", -1, &stmt, NULL);
        sqlite3_bind_int(stmt, 1, object_id);
        rc = sqlite3_step(stmt);
        check_rc(rc);

        // if name (or any other column) is NULL, this object_id is not in the database
        if (sqlite3_column_type(stmt, 0) == SQLITE_NULL)
            throw MovingTargetError("object_id " + std::to_string(target.object_id()) + " not found");

        // otherwise, loop through the names and add them to our object
        while (rc == SQLITE_ROW)
        {
            target.add_name(
                string(reinterpret_cast<const char *>(sqlite3_column_text(stmt, 0))),
                (bool)sqlite3_column_int(stmt, 1));
            rc = sqlite3_step(stmt);
            check_rc(rc);
        }

        sqlite3_finalize(stmt);

        target.object_id(object_id);
        return target;
    }

    MovingTarget SBSearchDatabaseSqlite3::get_moving_target(const string &name)
    {
        error_if_closed();
        sqlite3_stmt *stmt;
        int rc;

        sqlite3_prepare_v2(db, "SELECT object_id FROM moving_targets WHERE name=? LIMIT 1", -1, &stmt, NULL);
        sqlite3_bind_text(stmt, 1, name.c_str(), -1, SQLITE_TRANSIENT);
        rc = sqlite3_step(stmt);
        check_rc(rc);

        // if object_id is NULL, this name is not in the database
        if (sqlite3_column_type(stmt, 0) == SQLITE_NULL)
            throw MovingTargetError("name " + name + " not found");

        // otherwise, return the target based on object_id
        int object_id = sqlite3_column_int(stmt, 0);
        sqlite3_finalize(stmt);
        return get_moving_target(object_id);
    }

    void SBSearchDatabaseSqlite3::add_observation(Observation &observation)
    {
        error_if_closed();

        int rc;
        int64 observation_id;
        sqlite3_stmt *statement;

        if (observation.terms().size() == 0)
            throw std::runtime_error("Observation is missing index terms.");

        if (observation.observation_id() == UNDEFINED_OBSID)
        {
            // insert row and update observation object with observation_id
            rc = sqlite3_prepare_v2(db, "INSERT INTO observations VALUES (NULL, ?, ?, ?, ?, ?, ?) RETURNING observation_id;", -1, &statement, NULL);
            check_rc(rc);
        }
        else
        {
            // update existing observation
            rc = sqlite3_prepare_v2(db, "UPDATE observations SET source=?, product_id=?, mjd_start=?, mjd_stop=?, fov=?, terms=? WHERE observation_id=?;", -1, &statement, NULL);
            rc += sqlite3_bind_int64(statement, 7, observation.observation_id());
        }

        rc += sqlite3_bind_text(statement, 1, observation.source().c_str(), -1, SQLITE_TRANSIENT);
        rc += sqlite3_bind_text(statement, 2, observation.product_id().c_str(), -1, SQLITE_TRANSIENT);
        rc += sqlite3_bind_double(statement, 3, observation.mjd_start());
        rc += sqlite3_bind_double(statement, 4, observation.mjd_stop());
        rc += sqlite3_bind_text(statement, 5, observation.fov().c_str(), -1, SQLITE_TRANSIENT);
        rc += sqlite3_bind_text(statement, 6, observation.terms().c_str(), -1, SQLITE_TRANSIENT);

        rc = sqlite3_step(statement);

        if (observation.observation_id() == UNDEFINED_OBSID)
            observation.observation_id(sqlite3_column_int64(statement, 0));
        else if (rc != SQLITE_DONE)
        {
            Logger::error() << sqlite3_errmsg(db) << endl;
            throw std::runtime_error("Error updating observation in database");
        }

        sqlite3_finalize(statement);
    }

    Observation SBSearchDatabaseSqlite3::get_observation(const int64 observation_id)
    {
        error_if_closed();

        sqlite3_stmt *statement;

        sqlite3_prepare_v2(db, "SELECT source, product_id, mjd_start, mjd_stop, fov, terms FROM observations WHERE observation_id = ?;", -1, &statement, NULL);
        sqlite3_bind_int64(statement, 1, observation_id);

        int rc = sqlite3_step(statement);
        check_rc(rc);

        if (rc != SQLITE_ROW)
            throw std::runtime_error("No matching observation.");

        string source((char *)sqlite3_column_text(statement, 0));
        string product_id((char *)sqlite3_column_text(statement, 1));
        double mjd_start = sqlite3_column_double(statement, 2);
        double mjd_stop = sqlite3_column_double(statement, 3);
        string fov((char *)sqlite3_column_text(statement, 4));
        string terms((char *)sqlite3_column_text(statement, 5));

        sqlite3_finalize(statement);

        return Observation(source, product_id, mjd_start, mjd_stop, fov, terms, observation_id);
    }

    vector<Observation> SBSearchDatabaseSqlite3::find_observations(vector<string> query_terms)
    {
        // query_terms may be spatial-temporal, just spatial, or just temporal.
        error_if_closed();

        char *error_message;
        char statement[MAXIMUM_QUERY_CLAUSE_LENGTH + 100];
        char *statement_end = statement;
        int count = 0;
        string term_string;
        std::set<int64> approximate_matches;

        auto collect_found_rowids = [](void *found_ptr, int count, char **data, char **columns)
        {
            static int total = 0;
            total += count;

            std::set<int64> *found = static_cast<std::set<int64> *>(found_ptr);

            for (int i = 0; i < count; i++)
                found->insert(std::strtoll(data[i], NULL, 10));

            return 0;
        };

        // Query database with terms, but not too many at once
        statement_end = stpcpy(statement, "SELECT rowid FROM observations_terms_index WHERE terms MATCH '"); // if this string's length changes, edit lines marked ** below
        auto term = query_terms.begin();
        while (term != query_terms.end())
        {
            if (statement_end != (statement + 62)) // ** match base statement length
            {
                // this is not the first term in the list: append OR
                statement_end = stpcpy(statement_end, " OR ");
            }
            // append the term using quotes
            statement_end = stpcpy(statement_end, "\"");
            statement_end = stpcpy(statement_end, (*term).c_str());
            statement_end = stpcpy(statement_end, "\"");
            term++;

            if (((statement_end - statement) > MAXIMUM_QUERY_CLAUSE_LENGTH) | (term == query_terms.end()))
            {
                // we have enough terms or have exhausted them all
                strcpy(statement_end, "';");
                sqlite3_exec(db, statement, collect_found_rowids, &approximate_matches, &error_message);
                check_sql(error_message);
                statement_end = statement + 62; // ** match base statement length
            }

            count++;
        }

        Logger::debug() << "Searched " << count << " of " << query_terms.size() << " query terms."
                        << endl;

        return get_observations(approximate_matches.begin(), approximate_matches.end());
    }

    void SBSearchDatabaseSqlite3::check_rc(const int rc)
    {
        if ((rc != SQLITE_OK) & (rc != SQLITE_ROW) & (rc != SQLITE_DONE))
        {
            Logger::error() << "sqlite3 error (" << rc << ") " << sqlite3_errmsg(db) << endl;
            throw std::runtime_error("sqlite3 error");
        }
    }

    void SBSearchDatabaseSqlite3::check_sql(char *error_message)
    {
        error_if_closed();

        if (error_message != NULL)
        {
            Logger::error() << error_message << endl;
            sqlite3_free(error_message);
            throw std::runtime_error("\nSQL error\n");
        }
    }

    void SBSearchDatabaseSqlite3::error_if_closed()
    {
        if (db == NULL)
            throw std::runtime_error("Database is not open.");
    }
}