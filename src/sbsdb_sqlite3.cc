#include "config.h"

#include <iostream>
#include <stdexcept>
#include <string>

#include <sqlite3.h>

#include "logging.h"
#include "observation.h"
#include "sbsdb.h"
#include "sbsdb_sqlite3.h"

using std::endl;

namespace sbsearch
{
    SBSearchDatabaseSqlite3::SBSearchDatabaseSqlite3(const char *filename)
    {
        int rc = sqlite3_open(filename, &db);
        if (rc != SQLITE_OK)
        {
            Logger::error() << "Error opening database: " << sqlite3_errmsg(db) << endl;
            throw std::runtime_error("Error opening database");
        }
        else
            Logger::info() << "Opened sqlite3 database " << filename << endl;
        execute_sql("PRAGMA temp_store_directory = './';");
    }

    SBSearchDatabaseSqlite3::~SBSearchDatabaseSqlite3()
    {
        close();
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

        create_observations_indices();

        // configuration table and defaults
        execute_sql(R"(
CREATE TABLE IF NOT EXISTS configuration (
  parameter TEXT NOT NULL UNIQUE,
  value TEXT NOT NULL
);
INSERT OR IGNORE INTO configuration VALUES ('max_spatial_cells', '8');
INSERT OR IGNORE INTO configuration VALUES ('max_spatial_level', '4');
INSERT OR IGNORE INTO configuration VALUES ('min_spatial_level', '12');
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
        sqlite3_exec(db, statement, callback, callback_arg, &error_message);
        check_sql(error_message);
    }

    double SBSearchDatabaseSqlite3::get_double(const char *statement)
    {
        double value;
        execute_sql(statement, get_single_value_callback<double>, &value);
        return value;
    }

    int SBSearchDatabaseSqlite3::get_int(const char *statement)
    {
        int value;
        execute_sql(statement, get_single_value_callback<int>, &value);
        return value;
    }

    int64 SBSearchDatabaseSqlite3::get_int64(const char *statement)
    {
        int64 value;
        execute_sql(statement, get_single_value_callback<int64>, &value);
        return value;
    }

    string SBSearchDatabaseSqlite3::get_string(const char *statement)
    {
        char *value;
        execute_sql(statement, get_single_value_callback<char *>, value);
        return string(value);
    }

    void SBSearchDatabaseSqlite3::indexer_options(Indexer::Options options)
    {
        error_if_closed();

        char *error_message = NULL;
        int rc;
        sqlite3_stmt *statement;

        std::vector<std::string> parameters = {"max_spatial_cells",
                                               "max_spatial_level",
                                               "min_spatial_level",
                                               "temporal_resolution"};
        std::vector<std::string> values = {std::to_string(options.max_spatial_cells()),
                                           std::to_string(options.max_spatial_level()),
                                           std::to_string(options.min_spatial_level()),
                                           std::to_string(options.temporal_resolution())};
        for (int i = 0; i < parameters.size(); i++)
        {
            sqlite3_prepare_v2(db, "UPDATE configuration SET parameter=?1, value=?2 WHERE parameter=?1;", -1, &statement, NULL);
            sqlite3_bind_text(statement, 1, parameters[i].c_str(), parameters[i].size(), SQLITE_STATIC);
            sqlite3_bind_text(statement, 2, values[i].c_str(), values[i].size(), SQLITE_STATIC);
            rc = sqlite3_step(statement);
            check_sql(error_message);
            sqlite3_finalize(statement);
        }
    }

    std::pair<double, double> SBSearchDatabaseSqlite3::date_range(string source)
    {
        double mjd_start, mjd_stop;
        if (source == "")
        {
            mjd_start = get_double("SELECT MIN(mjd_start) FROM observations;");
            mjd_stop = get_double("SELECT MAX(mjd_stop) FROM observations;");
        }
        else
        {
            char *error_message = NULL;
            sqlite3_stmt *statement;
            int rc;
            rc = sqlite3_prepare_v2(db, "SELECT MIN(mjd_start), MAX(mjd_stop) FROM observations WHERE source=?;", -1, &statement, NULL);
            check_rc(rc);
            rc += sqlite3_bind_text(statement, 1, source.c_str(), -1, SQLITE_TRANSIENT);
            check_sql(error_message);

            rc = sqlite3_step(statement);
            check_sql(error_message);

            mjd_start = sqlite3_column_double(statement, 0);
            mjd_stop = sqlite3_column_double(statement, 1);

            sqlite3_finalize(statement);
        }

        return std::pair<double, double>(mjd_start, mjd_stop);
    }

    void SBSearchDatabaseSqlite3::add_observation(Observation &observation)
    {
        error_if_closed();

        char *error_message = NULL;
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
        check_sql(error_message);

        rc = sqlite3_step(statement);
        check_sql(error_message);

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

        char *error_message = 0;
        sqlite3_stmt *statement;

        sqlite3_prepare_v2(db, "SELECT source, product_id, mjd_start, mjd_stop, fov, terms FROM observations WHERE observation_id = ?;", -1, &statement, NULL);
        sqlite3_bind_int64(statement, 1, observation_id);
        check_sql(error_message);

        int rc = sqlite3_step(statement);
        check_sql(error_message);

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