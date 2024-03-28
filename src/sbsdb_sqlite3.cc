#include "config.h"

#include <algorithm>
#include <ctime>
#include <iostream>
#include <stdexcept>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <sqlite3.h>

#include "ephemeris.h"
#include "found.h"
#include "exceptions.h"
#include "logging.h"
#include "moving_target.h"
#include "observation.h"
#include "observatory.h"
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
        int flags = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE;
        if (filename == ":memory:")
            flags |= SQLITE_OPEN_MEMORY;

        int rc = sqlite3_open_v2(filename.c_str(), &db, flags, NULL);
        if (rc != SQLITE_OK)
        {
            Logger::error() << "Error opening database: " << sqlite3_errmsg(db) << endl;
            throw std::runtime_error("Error opening database");
        }
        else
            Logger::info() << "Opened sqlite3 database " << filename << endl;
        execute_sql("PRAGMA temp_store_directory = './';");

        setup_tables();
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
  observatory TEXT NOT NULL,
  product_id TEXT NOT NULL,
  mjd_start FLOAT NOT NULL,
  mjd_stop FLOAT NOT NULL,
  fov TEXT NOT NULL,
  terms TEXT NOT NULL
);
)");
        create_observations_indices();

        // The moving_target_id defines a unique object.
        execute_sql(R"(
CREATE TABLE IF NOT EXISTS moving_targets (
  moving_targets_row_id INTEGER PRIMARY KEY,
  moving_target_id INTEGER NOT NULL,
  name TEXT NOT NULL,
  small_body BOOLEAN NOT NULL,
  primary_id BOOLEAN NOT NULL
);
CREATE UNIQUE INDEX IF NOT EXISTS idx_moving_target_primary_id ON moving_targets(moving_target_id, primary_id) WHERE primary_id=TRUE;
CREATE UNIQUE INDEX IF NOT EXISTS idx_moving_target_name_small_body ON moving_targets(name, small_body);
CREATE INDEX IF NOT EXISTS idx_moving_target_moving_target_id ON moving_targets(moving_target_id);
)");

        execute_sql(R"(
CREATE TABLE IF NOT EXISTS observatories (
  name TEXT UNIQUE NOT NULL,
  longitude FLOAT NOT NULL,
  rho_cos_phi FLOAT NOT NULL,
  rho_sin_phi FLOAT NOT NULL
);
)");

        execute_sql(R"(
CREATE TABLE IF NOT EXISTS ephemerides (
  ephemeris_id INTEGER PRIMARY KEY,
  moving_target_id INTEGER NOT NULL,
  mjd FLOAT NOT NULL,
  tmtp FLOAT NOT NULL,
  ra FLOAT NOT NULL,
  dec FLOAT NOT NULL,
  unc_a FLOAT NOT NULL,
  unc_b FLOAT NOT NULL,
  unc_theta FLOAT NOT NULL,
  rh FLOAT NOT NULL,
  delta FLOAT NOT NULL,
  phase FLOAT NOT NULL,
  selong FLOAT NOT NULL,
  true_anomaly FLOAT NOT NULL,
  sangle FLOAT NOT NULL,
  vangle FLOAT NOT NULL,
  vmag FLOAT NOT NULL,
  retrieved TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_ephemerides_moving_target_id ON ephemerides(moving_target_id);
)");

        execute_sql(R"(
CREATE TABLE IF NOT EXISTS found (
  found_id INTEGER PRIMARY KEY,
  observation_id INTEGER NOT NULL,
  moving_target_id INTEGER NOT NULL,
  mjd FLOAT NOT NULL,
  tmtp FLOAT NOT NULL,
  ra FLOAT NOT NULL,
  dec FLOAT NOT NULL,
  unc_a FLOAT NOT NULL,
  unc_b FLOAT NOT NULL,
  unc_theta FLOAT NOT NULL,
  rh FLOAT NOT NULL,
  delta FLOAT NOT NULL,
  phase FLOAT NOT NULL,
  selong FLOAT NOT NULL,
  true_anomaly FLOAT NOT NULL,
  sangle FLOAT NOT NULL,
  vangle FLOAT NOT NULL,
  vmag FLOAT NOT NULL,
  saved TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_found_observation_id ON found(observation_id);
CREATE UNIQUE INDEX IF NOT EXISTS idx_found_moving_target_id ON found(moving_target_id, observation_id);
)");

        // configuration table and defaults
        execute_sql(R"(
CREATE TABLE IF NOT EXISTS configuration (
  parameter TEXT NOT NULL UNIQUE,
  value TEXT NOT NULL
);
INSERT OR IGNORE INTO configuration VALUES ('max_spatial_index_cells', '8');
INSERT OR IGNORE INTO configuration VALUES ('max_spatial_level', '12');
INSERT OR IGNORE INTO configuration VALUES ('min_spatial_level', '4');
INSERT OR IGNORE INTO configuration VALUES ('temporal_resolution', '1');
INSERT OR IGNORE INTO configuration VALUES ('database version', ')" SBSEARCH_DATABASE_VERSION "');");

        Logger::debug() << "Database tables are set." << endl;
    }

    void SBSearchDatabaseSqlite3::create_observations_indices()
    {
        Logger::info() << "Checking for observations indices." << std::endl;
        // Consider a tokenizer that includes hyphens.

        error_if_closed();

        char *error_message = NULL;
        int rc = sqlite3_exec(db, R"(
  CREATE VIRTUAL TABLE observations_terms_index
    USING fts5(terms, content='observations', content_rowid='observation_id');
)",
                              NULL, NULL, &error_message);

        // if the create table statement errors, do nothing, otherwise report
        // the error
        if ((rc != SQLITE_OK) & (rc != SQLITE_ROW) & (rc != SQLITE_DONE))
        {
            if (string(error_message) != "table observations_terms_index already exists")
                check_sql(error_message);
            sqlite3_free(error_message);
        }
        else
        {
            // the table was created, add any terms that already exist in the database
            execute_sql("INSERT INTO observations_terms_index (rowid, terms)"
                        " SELECT observation_id, terms FROM observations;");
            Logger::info() << "Created observation term index." << std::endl;
        }

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

        execute_sql(
            "CREATE INDEX IF NOT EXISTS idx_observations_mjd_start ON observations(mjd_start);\n"
            "CREATE INDEX IF NOT EXISTS idx_observations_mjd_stop ON observations(mjd_stop);\n"
            "CREATE INDEX IF NOT EXISTS idx_observations_source_mjd_start ON observations(source, mjd_start);\n"
            "CREATE INDEX IF NOT EXISTS idx_observations_source_mjd_stop ON observations(source, mjd_stop);\n"
            "CREATE UNIQUE INDEX IF NOT EXISTS idx_observations_source_product_id ON observations(source, product_id);\n"
            "CREATE INDEX IF NOT EXISTS idx_observations_product_id ON observations(product_id);\n");
    }

    void SBSearchDatabaseSqlite3::drop_observations_indices()
    {
        Logger::info() << "Dropping observations indices." << std::endl;
        execute_sql(
            "DROP TRIGGER IF EXISTS observations_insert;"
            "DROP TRIGGER IF EXISTS observations_delete;"
            "DROP TRIGGER IF EXISTS observations_update;"
            "DROP TABLE IF EXISTS observations_terms_index;"
            "DROP INDEX IF EXISTS idx_observations_mjd_start;"
            "DROP INDEX IF EXISTS idx_observations_mjd_stop;");
        Logger::info() << "Observations indices dropped." << std::endl;
    }

    void SBSearchDatabaseSqlite3::execute_sql(const char *statement) const
    {
        execute_sql(statement, NULL, NULL);
    }

    void SBSearchDatabaseSqlite3::execute_sql(const char *statement, int (*callback)(void *, int, char **, char **), void *callback_arg) const
    {
        error_if_closed();

        char *error_message = NULL;
        int rc = sqlite3_exec(db, statement, callback, callback_arg, &error_message);
        check_sql(error_message);
    }

    double *SBSearchDatabaseSqlite3::get_double(const char *statement) const
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

    int *SBSearchDatabaseSqlite3::get_int(const char *statement) const
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

    int64 *SBSearchDatabaseSqlite3::get_int64(const char *statement) const
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

    string *SBSearchDatabaseSqlite3::get_string(const char *statement) const
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

        vector<string> parameters = {"max_spatial_index_cells",
                                     "max_spatial_level",
                                     "min_spatial_level",
                                     "temporal_resolutionstd::"};
        vector<string> values = {std::to_string(options.max_spatial_index_cells()),
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

    std::pair<double *, double *> SBSearchDatabaseSqlite3::observation_date_range(const string &source) const
    {
        double *mjd_start = new double;
        double *mjd_stop = new double;
        if (source.empty() | (source == ""))
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

    void SBSearchDatabaseSqlite3::add_moving_target(MovingTarget &target) const
    {
        error_if_closed();

        int moving_target_id = target.moving_target_id();
        if (moving_target_id == UNDEF_MOVING_TARGET_ID)
        {
            // new objects use db largest moving_target_id + 1, or 1 if there are no objects
            moving_target_id = *get_int("SELECT IFNULL(MAX(moving_target_id), 0) + 1 FROM moving_targets");
            target.moving_target_id(moving_target_id);
        }

        int rc;
        sqlite3_stmt *stmt;
        sqlite3_prepare_v2(db, "SELECT COUNT() FROM moving_targets WHERE moving_target_id=?", -1, &stmt, NULL);
        sqlite3_bind_int(stmt, 1, target.moving_target_id());
        rc = sqlite3_step(stmt);
        check_rc(rc);
        int count = sqlite3_column_int(stmt, 0);
        if (count != 0)
            throw MovingTargetError("moving target id " +
                                    std::to_string(target.moving_target_id()) +
                                    " already exists");
        sqlite3_finalize(stmt);

        Logger::info() << "Add moving target " << target << endl;
        try
        {
            execute_sql("BEGIN TRANSACTION;");
            add_moving_target_name(moving_target_id, target.designation(), target.small_body(), true);
            for (const string &name : target.alternate_names())
                add_moving_target_name(target.moving_target_id(), name, target.small_body(), false);
            execute_sql("END TRANSACTION;");
        }
        catch (const MovingTargetError &err)
        {
            execute_sql("ROLLBACK TRANSACTION;");
            throw err;
        }
        Logger::info() << target << " added to database." << std::endl;
    }

    void SBSearchDatabaseSqlite3::add_moving_target_name(const int moving_target_id,
                                                         const string &name,
                                                         const bool small_body,
                                                         const bool primary_id) const
    {
        Logger::debug() << "Add moving target " << name
                        << " (ID=" << moving_target_id
                        << "; small body=" << (small_body ? "true" : "false")
                        << "; primary=" << (primary_id ? "true" : "false")
                        << ")" << std::endl;

        error_if_closed();
        int rc;
        sqlite3_stmt *stmt;
        sqlite3_prepare_v2(db, "INSERT INTO moving_targets (moving_target_id, name, small_body, primary_id) VALUES (?, ?, ?, ?)", -1, &stmt, NULL);
        sqlite3_bind_int(stmt, 1, moving_target_id);
        sqlite3_bind_text(stmt, 2, name.c_str(), -1, SQLITE_TRANSIENT);
        sqlite3_bind_int(stmt, 3, small_body);
        sqlite3_bind_int(stmt, 4, primary_id);
        rc = sqlite3_step(stmt);
        if (rc == SQLITE_CONSTRAINT)
        {
            const char *message = sqlite3_errmsg(db);
            throw MovingTargetError("target uniqueness violated: " + string(message));
        }
        check_rc(rc);
        sqlite3_finalize(stmt);
    };

    void SBSearchDatabaseSqlite3::remove_moving_target(const MovingTarget &target) const
    {
        error_if_closed();
        int rc;
        sqlite3_stmt *stmt;

        try
        {
            execute_sql("BEGIN TRANSACTION;");
            sqlite3_prepare_v2(db, "SELECT COUNT() FROM moving_targets WHERE moving_target_id=? AND small_body=?", -1, &stmt, NULL);
            sqlite3_bind_int(stmt, 1, target.moving_target_id());
            sqlite3_bind_int(stmt, 2, target.small_body());
            rc = sqlite3_step(stmt);
            check_rc(rc);
            int count = sqlite3_column_int(stmt, 0);
            if (count == 0)
                throw MovingTargetError("moving target id " +
                                        std::to_string(target.moving_target_id()) +
                                        " with small body flag " +
                                        (target.small_body() ? "true" : "false") +
                                        " not found");
            sqlite3_finalize(stmt);

            Logger::info() << "Remove moving target " << target << endl;
            sqlite3_prepare_v2(db, "DELETE FROM moving_targets WHERE moving_target_id=? AND small_body=?", -1, &stmt, NULL);
            sqlite3_bind_int(stmt, 1, target.moving_target_id());
            sqlite3_bind_int(stmt, 2, target.small_body());
            rc = sqlite3_step(stmt);
            check_rc(rc);
            sqlite3_finalize(stmt);
            execute_sql("END TRANSACTION;");
        }
        catch (const MovingTargetError &err)
        {
            execute_sql("ROLLBACK TRANSACTION;");
            throw err;
        }
        Logger::info() << target << " removed from database." << std::endl;
    }

    void SBSearchDatabaseSqlite3::update_moving_target(const MovingTarget &target) const
    {
        Logger::info() << "Update moving target " << target << endl;
        remove_moving_target(target);
        MovingTarget copy(target);
        add_moving_target(copy);
        Logger::info() << target << " updated." << std::endl;
    }

    MovingTarget SBSearchDatabaseSqlite3::get_moving_target(const int moving_target_id) const
    {
        error_if_closed();
        MovingTarget target;

        sqlite3_stmt *stmt;
        int rc;

        sqlite3_prepare_v2(db, "SELECT name, small_body, primary_id FROM moving_targets WHERE moving_target_id=?", -1, &stmt, NULL);
        sqlite3_bind_int(stmt, 1, moving_target_id);
        rc = sqlite3_step(stmt);
        check_rc(rc);

        // if name (or any other column) is NULL, this moving_target_id is not in the database
        if (sqlite3_column_type(stmt, 0) == SQLITE_NULL)
            throw MovingTargetError("moving target id " +
                                    std::to_string(target.moving_target_id()) +
                                    " not found");

        target.small_body((bool)sqlite3_column_int(stmt, 1));

        // otherwise, loop through the names and add them to our object
        while (rc == SQLITE_ROW)
        {
            string name(reinterpret_cast<const char *>(sqlite3_column_text(stmt, 0)));
            bool primary = (bool)sqlite3_column_int(stmt, 2);

            if (primary)
                target.designation(name);
            else
                target.add_name(name);

            rc = sqlite3_step(stmt);
            check_rc(rc);
        }

        sqlite3_finalize(stmt);

        target.moving_target_id(moving_target_id);
        return target;
    }

    MovingTarget SBSearchDatabaseSqlite3::get_moving_target(const string &name, const bool small_body) const
    {
        error_if_closed();
        sqlite3_stmt *stmt;
        int rc;

        sqlite3_prepare_v2(db, "SELECT moving_target_id FROM moving_targets WHERE name=? AND small_body=? LIMIT 1", -1, &stmt, NULL);
        sqlite3_bind_text(stmt, 1, name.c_str(), -1, SQLITE_TRANSIENT);
        sqlite3_bind_int(stmt, 2, small_body);
        rc = sqlite3_step(stmt);
        check_rc(rc);

        // if moving_target_id is NULL, this name-small body combo is not in the
        // database
        if (sqlite3_column_type(stmt, 0) == SQLITE_NULL)
            return MovingTarget(name, small_body);

        // otherwise, return the target based on moving_target_id
        int moving_target_id = sqlite3_column_int(stmt, 0);
        sqlite3_finalize(stmt);
        return get_moving_target(moving_target_id);
    }

    vector<MovingTarget> SBSearchDatabaseSqlite3::get_all_moving_targets() const
    {
        error_if_closed();
        sqlite3_stmt *stmt;
        int rc;

        sqlite3_prepare_v2(db, "SELECT DISTINCT(moving_target_id) FROM moving_targets", -1, &stmt, NULL);
        rc = sqlite3_step(stmt);
        check_rc(rc);

        vector<MovingTarget> targets;
        while (rc == SQLITE_ROW)
        {
            targets.push_back(get_moving_target(sqlite3_column_int(stmt, 0)));
            rc = sqlite3_step(stmt);
            check_rc(rc);
        }
        sqlite3_finalize(stmt);
        return targets;
    }

    void SBSearchDatabaseSqlite3::add_observatory(const string &name, const Observatory &observatory) const
    {
        error_if_closed();

        // do not add anything if this name is already in the database
        Observatory in_database;
        try
        {
            in_database = get_observatory(name);
        }
        catch (const ObservatoryError &)
        {
            Logger::info() << "Adding observatory for " << name << " to the database." << std::endl;
        }
        if (in_database != Observatory())
            throw ObservatoryError(name + " already in the database");

        int rc;
        sqlite3_stmt *stmt;
        sqlite3_prepare_v2(db, R"(
INSERT INTO observatories (
  name, longitude, rho_cos_phi, rho_sin_phi
) VALUES (
  ?, ?, ?, ?
);)",
                           -1, &stmt, NULL);
        sqlite3_bind_text(stmt, 1, name.c_str(), -1, SQLITE_TRANSIENT);
        sqlite3_bind_double(stmt, 2, observatory.longitude);
        sqlite3_bind_double(stmt, 3, observatory.rho_cos_phi);
        sqlite3_bind_double(stmt, 4, observatory.rho_sin_phi);
        rc = sqlite3_step(stmt);
        check_rc(rc);
        sqlite3_finalize(stmt);
    }

    const Observatory SBSearchDatabaseSqlite3::get_observatory(const string &name) const
    {
        int rc;
        sqlite3_stmt *stmt;

        sqlite3_prepare_v2(db, R"(
SELECT longitude, rho_cos_phi, rho_sin_phi
FROM observatories
WHERE name = ?;
)",
                           -1, &stmt, NULL);
        sqlite3_bind_text(stmt, 1, name.c_str(), -1, SQLITE_TRANSIENT);
        rc = sqlite3_step(stmt);
        check_rc(rc);

        if (sqlite3_column_type(stmt, 0) == SQLITE_NULL)
            throw ObservatoryError(name + " not found");

        Observatory observatory{sqlite3_column_double(stmt, 0),
                                sqlite3_column_double(stmt, 1),
                                sqlite3_column_double(stmt, 2)};
        sqlite3_finalize(stmt);
        return observatory;
    }

    const Observatories SBSearchDatabaseSqlite3::get_observatories() const
    {
        int rc;
        sqlite3_stmt *stmt;
        Observatories observatories;

        sqlite3_prepare_v2(db, "SELECT name, longitude, rho_cos_phi, rho_sin_phi FROM observatories",
                           -1, &stmt, NULL);
        rc = sqlite3_step(stmt);
        check_rc(rc);

        while (rc == SQLITE_ROW)
        {
            const string name((char *)sqlite3_column_text(stmt, 0));
            const Observatory observatory{sqlite3_column_double(stmt, 1),
                                          sqlite3_column_double(stmt, 2),
                                          sqlite3_column_double(stmt, 3)};
            observatories[name] = observatory;

            rc = sqlite3_step(stmt);
            check_rc(rc);
        }

        sqlite3_finalize(stmt);

        return observatories;
    }

    void SBSearchDatabaseSqlite3::remove_observatory(const string &name) const
    {
        error_if_closed();
        int rc;
        sqlite3_stmt *stmt;

        Logger::info() << "Removing observatory for name " << name << endl;
        sqlite3_prepare_v2(db, "DELETE FROM observatories WHERE name=?", -1, &stmt, NULL);
        sqlite3_bind_text(stmt, 1, name.c_str(), -1, SQLITE_TRANSIENT);
        rc = sqlite3_step(stmt);
        check_rc(rc);
        sqlite3_finalize(stmt);
    }

    const vector<string> SBSearchDatabaseSqlite3::get_sources() const
    {
        error_if_closed();

        sqlite3_stmt *statement;

        sqlite3_prepare_v2(db, "SELECT DISTINCT(source) FROM observations", -1, &statement, NULL);

        int rc = sqlite3_step(statement);
        check_rc(rc);
        vector<string> sources;
        while (rc == SQLITE_ROW)
        {
            sources.push_back((char *)sqlite3_column_text(statement, 0));
            rc = sqlite3_step(statement);
            check_rc(rc);
        }
        sqlite3_finalize(statement);
        return sources;
    }

    void SBSearchDatabaseSqlite3::add_ephemeris(Ephemeris &eph) const
    {
        error_if_closed();

        // verify that the moving target ID exists in the database
        MovingTarget target = get_moving_target(eph.target().moving_target_id()); // throws MovingTargetError if not found
        if (target != eph.target())
            throw MovingTargetError("Ephemeris target does not match database copy");

        Logger::info()
            << "Adding " << std::to_string(eph.num_vertices()) << " ephemeris epochs for target " << eph.target().designation() << " (moving_target_id=" << eph.target().moving_target_id() << ")." << endl;

        char now[32];
        std::time_t time_now = std::time(nullptr);
        std::strftime(now, 32, "%F %T", std::gmtime(&time_now));

        int rc;
        sqlite3_stmt *stmt;
        execute_sql("BEGIN TRANSACTION;");
        for (const Ephemeris::Datum row : eph.data())
        {
            sqlite3_prepare_v2(db, R"(
INSERT INTO ephemerides (
  moving_target_id, mjd, tmtp,
  ra, dec, unc_a, unc_b, unc_theta,
  rh, delta, phase, selong, true_anomaly,
  sangle, vangle, vmag, retrieved
) VALUES (
  ?, ?, ?, ?, ?, ?, ?, ?,
  ?, ?, ?, ?, ?, ?, ?, ?, ?
);)",
                               -1, &stmt, NULL);
            sqlite3_bind_int(stmt, 1, eph.target().moving_target_id());
            sqlite3_bind_double(stmt, 2, row.mjd);
            sqlite3_bind_double(stmt, 3, row.tmtp);
            sqlite3_bind_double(stmt, 4, row.ra);
            sqlite3_bind_double(stmt, 5, row.dec);
            sqlite3_bind_double(stmt, 6, row.unc_a);
            sqlite3_bind_double(stmt, 7, row.unc_b);
            sqlite3_bind_double(stmt, 8, row.unc_theta);
            sqlite3_bind_double(stmt, 9, row.rh);
            sqlite3_bind_double(stmt, 10, row.delta);
            sqlite3_bind_double(stmt, 11, row.phase);
            sqlite3_bind_double(stmt, 12, row.selong);
            sqlite3_bind_double(stmt, 13, row.true_anomaly);
            sqlite3_bind_double(stmt, 14, row.sangle);
            sqlite3_bind_double(stmt, 15, row.vangle);
            sqlite3_bind_double(stmt, 16, row.vmag);
            sqlite3_bind_text(stmt, 17, now, -1, SQLITE_STATIC);
            rc = sqlite3_step(stmt);
            check_rc(rc);
            sqlite3_finalize(stmt);
        }
        execute_sql("END TRANSACTION;");
    }

    Ephemeris SBSearchDatabaseSqlite3::get_ephemeris(const MovingTarget target, double mjd_start, double mjd_stop) const
    {
        int rc;
        sqlite3_stmt *stmt;
        Ephemeris::Data data;

        sqlite3_prepare_v2(db, "SELECT COUNT() FROM ephemerides WHERE moving_target_id=? AND mjd >= ? and mjd <= ?;", -1, &stmt, NULL);
        sqlite3_bind_int(stmt, 1, target.moving_target_id());
        sqlite3_bind_double(stmt, 2, mjd_start);
        sqlite3_bind_double(stmt, 3, mjd_stop);
        rc = sqlite3_step(stmt);
        check_rc(rc);
        int count = sqlite3_column_int(stmt, 0);
        sqlite3_finalize(stmt);
        Logger::debug() << "Reading " << count << " ephemeris epochs from database for " << target.designation() << " (moving_target_id=" << target.moving_target_id() << ")." << endl;

        data.reserve(count);
        sqlite3_prepare_v2(db, R"(
SELECT
    mjd, tmtp, ra, dec, unc_a, unc_b, unc_theta,
    rh, delta, phase, selong, true_anomaly,
    sangle, vangle, vmag
FROM ephemerides
WHERE moving_target_id=? AND mjd >= ? and mjd <= ?;)",
                           -1, &stmt, NULL);
        sqlite3_bind_int(stmt, 1, target.moving_target_id());
        sqlite3_bind_double(stmt, 2, mjd_start);
        sqlite3_bind_double(stmt, 3, mjd_stop);

        for (int i = 0; i < count; i++)
        {
            Ephemeris::Datum d;
            rc = sqlite3_step(stmt);
            check_rc(rc);
            d.mjd = sqlite3_column_double(stmt, 0);
            d.tmtp = sqlite3_column_double(stmt, 1);
            d.ra = sqlite3_column_double(stmt, 2);
            d.dec = sqlite3_column_double(stmt, 3);
            d.unc_a = sqlite3_column_double(stmt, 4);
            d.unc_b = sqlite3_column_double(stmt, 5);
            d.unc_theta = sqlite3_column_double(stmt, 6);
            d.rh = sqlite3_column_double(stmt, 7);
            d.delta = sqlite3_column_double(stmt, 8);
            d.phase = sqlite3_column_double(stmt, 9);
            d.selong = sqlite3_column_double(stmt, 10);
            d.true_anomaly = sqlite3_column_double(stmt, 11);
            d.sangle = sqlite3_column_double(stmt, 12);
            d.vangle = sqlite3_column_double(stmt, 13);
            d.vmag = sqlite3_column_double(stmt, 14);
            data.push_back(d);
        }
        sqlite3_finalize(stmt);

        return {target, data};
    }

    int SBSearchDatabaseSqlite3::remove_ephemeris(const MovingTarget target, double mjd_start, double mjd_stop) const
    {
        int rc;
        sqlite3_stmt *stmt;

        sqlite3_prepare_v2(db, "SELECT COUNT() FROM ephemerides WHERE moving_target_id=? AND mjd >= ? and mjd <= ?;", -1, &stmt, NULL);
        sqlite3_bind_int(stmt, 1, target.moving_target_id());
        sqlite3_bind_double(stmt, 2, mjd_start);
        sqlite3_bind_double(stmt, 3, mjd_stop);
        rc = sqlite3_step(stmt);
        check_rc(rc);
        int count = sqlite3_column_int(stmt, 0);
        sqlite3_finalize(stmt);
        Logger::info() << "Removing " << count << " ephemeris epochs from database for " << target.designation() << " (moving_target_id=" << target.moving_target_id() << ")." << endl;

        sqlite3_prepare_v2(db, "DELETE FROM ephemerides WHERE moving_target_id=? AND mjd >= ? and mjd <= ?;",
                           -1, &stmt, NULL);
        sqlite3_bind_int(stmt, 1, target.moving_target_id());
        sqlite3_bind_double(stmt, 2, mjd_start);
        sqlite3_bind_double(stmt, 3, mjd_stop);
        rc = sqlite3_step(stmt);
        check_rc(rc);
        sqlite3_finalize(stmt);

        return count;
    }

    std::pair<double *, double *> SBSearchDatabaseSqlite3::ephemeris_date_range() const
    {
        double *mjd_start = new double;
        double *mjd_stop = new double;

        mjd_start = get_double("SELECT MIN(mjd) FROM ephemerides;");
        mjd_stop = get_double("SELECT MAX(mjd) FROM ephemerides;");

        return std::pair<double *, double *>(std::move(mjd_start), std::move(mjd_stop));
    }

    void SBSearchDatabaseSqlite3::add_observation(Observation &observation) const
    {
        error_if_closed();

        int rc;
        int64 observation_id;
        sqlite3_stmt *statement;

        if (observation.terms().size() == 0)
            throw std::runtime_error("Observation is missing index terms.");

        int index = 0;
        if (observation.observation_id() == UNDEFINED_OBSID)
        {
            // insert row and update observation object with observation_id
            rc = sqlite3_prepare_v2(db, "INSERT INTO observations VALUES (NULL, ?, ?, ?, ?, ?, ?, ?) RETURNING observation_id;", -1, &statement, NULL);
        }
        else
        {
            // update existing observation
            sqlite3_prepare_v2(db, "INSERT OR REPLACE INTO observations VALUES (?, ?, ?, ?, ?, ?, ?, ?);", -1, &statement, NULL);
            sqlite3_bind_int64(statement, ++index, observation.observation_id());
        }

        sqlite3_bind_text(statement, ++index, observation.source().c_str(), -1, SQLITE_TRANSIENT);
        sqlite3_bind_text(statement, ++index, observation.observatory().c_str(), -1, SQLITE_TRANSIENT);
        sqlite3_bind_text(statement, ++index, observation.product_id().c_str(), -1, SQLITE_TRANSIENT);
        sqlite3_bind_double(statement, ++index, observation.mjd_start());
        sqlite3_bind_double(statement, ++index, observation.mjd_stop());
        sqlite3_bind_text(statement, ++index, observation.fov().c_str(), -1, SQLITE_TRANSIENT);
        sqlite3_bind_text(statement, ++index, observation.terms().c_str(), -1, SQLITE_TRANSIENT);

        rc = sqlite3_step(statement);
        check_rc(rc);

        if (observation.observation_id() == UNDEFINED_OBSID)
            observation.observation_id(sqlite3_column_int64(statement, 0));
        else if (rc != SQLITE_DONE)
        {
            Logger::error() << sqlite3_errmsg(db) << endl;
            throw std::runtime_error("Error updating observation in database");
        }

        sqlite3_finalize(statement);
    }

    Observation SBSearchDatabaseSqlite3::get_observation(const int64 observation_id) const
    {
        error_if_closed();

        sqlite3_stmt *statement;

        sqlite3_prepare_v2(db, "SELECT source, observatory, product_id, mjd_start, mjd_stop, fov, terms FROM observations WHERE observation_id = ?;", -1, &statement, NULL);
        sqlite3_bind_int64(statement, 1, observation_id);

        int rc = sqlite3_step(statement);
        check_rc(rc);

        if (rc != SQLITE_ROW)
            throw std::runtime_error("No matching observation.");

        string source((char *)sqlite3_column_text(statement, 0));
        string observatory((char *)sqlite3_column_text(statement, 1));
        string product_id((char *)sqlite3_column_text(statement, 2));
        double mjd_start = sqlite3_column_double(statement, 3);
        double mjd_stop = sqlite3_column_double(statement, 4);
        string fov((char *)sqlite3_column_text(statement, 5));
        string terms((char *)sqlite3_column_text(statement, 6));

        sqlite3_finalize(statement);

        return Observation(source, observatory, product_id, mjd_start, mjd_stop, fov, terms, observation_id);
    }

    void SBSearchDatabaseSqlite3::remove_observations(const double mjd_start, const double mjd_stop) const
    {
        error_if_closed();

        sqlite3_stmt *statement;

        sqlite3_prepare_v2(db, "DELETE FROM observations WHERE mjd_start >= ? AND mjd_stop <= ?", -1, &statement, NULL);
        sqlite3_bind_double(statement, 1, mjd_start);
        sqlite3_bind_double(statement, 2, mjd_stop);

        int rc = sqlite3_step(statement);
        check_rc(rc);
        sqlite3_finalize(statement);
    }

    void SBSearchDatabaseSqlite3::remove_observations(const string &source, const double mjd_start, const double mjd_stop) const
    {
        error_if_closed();

        sqlite3_stmt *statement;

        sqlite3_prepare_v2(db, "DELETE FROM observations WHERE source = ? AND mjd_start >= ? AND mjd_stop <= ?", -1, &statement, NULL);
        sqlite3_bind_text(statement, 1, source.c_str(), -1, SQLITE_TRANSIENT);
        sqlite3_bind_double(statement, 2, mjd_start);
        sqlite3_bind_double(statement, 3, mjd_stop);

        int rc = sqlite3_step(statement);
        check_rc(rc);
        sqlite3_finalize(statement);
    }

    int64 SBSearchDatabaseSqlite3::count_observations(const double mjd_start, const double mjd_stop) const
    {
        error_if_closed();

        sqlite3_stmt *statement;

        sqlite3_prepare_v2(db, "SELECT COUNT() FROM observations WHERE mjd_start >= ? AND mjd_stop <= ?;", -1, &statement, NULL);
        sqlite3_bind_double(statement, 1, mjd_start);
        sqlite3_bind_double(statement, 2, mjd_stop);

        int rc = sqlite3_step(statement);
        check_rc(rc);

        int64 count = sqlite3_column_int64(statement, 0);
        sqlite3_finalize(statement);

        return count;
    }

    int64 SBSearchDatabaseSqlite3::count_observations(const string &source, const double mjd_start, const double mjd_stop) const
    {
        if (source == "")
            return count_observations(mjd_start, mjd_stop);

        error_if_closed();

        sqlite3_stmt *statement;

        sqlite3_prepare_v2(db, "SELECT COUNT() FROM observations WHERE source = ? AND mjd_start >= ? AND mjd_stop <= ?;", -1, &statement, NULL);
        sqlite3_bind_text(statement, 1, source.c_str(), -1, SQLITE_TRANSIENT);
        sqlite3_bind_double(statement, 2, mjd_start);
        sqlite3_bind_double(statement, 3, mjd_stop);

        int rc = sqlite3_step(statement);
        check_rc(rc);

        int64 count = sqlite3_column_int64(statement, 0);
        sqlite3_finalize(statement);
        return count;
    }

    Observations SBSearchDatabaseSqlite3::find_observations(const double mjd_start, const double mjd_stop, const int64 limit, const int64 offset) const
    {
        error_if_closed();

        sqlite3_stmt *statement;

        sqlite3_prepare_v2(db, "SELECT observation_id, source, observatory, product_id, mjd_start, mjd_stop, fov, terms FROM observations WHERE mjd_start >= ? AND mjd_stop <= ? LIMIT ? OFFSET ?;", -1, &statement, NULL);
        sqlite3_bind_double(statement, 1, mjd_start);
        sqlite3_bind_double(statement, 2, mjd_stop);
        sqlite3_bind_int64(statement, 3, limit);
        sqlite3_bind_int64(statement, 4, offset);

        int rc = sqlite3_step(statement);
        check_rc(rc);

        Observations observations;
        observations.reserve(limit);
        while (rc == SQLITE_ROW)
        {
            observations.push_back({string((char *)sqlite3_column_text(statement, 1)),
                                    string((char *)sqlite3_column_text(statement, 2)),
                                    string((char *)sqlite3_column_text(statement, 3)),
                                    sqlite3_column_double(statement, 4),
                                    sqlite3_column_double(statement, 5),
                                    string((char *)sqlite3_column_text(statement, 6)),
                                    string((char *)sqlite3_column_text(statement, 7)),
                                    sqlite3_column_int64(statement, 0)});
            rc = sqlite3_step(statement);
            check_rc(rc);
        }

        sqlite3_finalize(statement);
        return observations;
    }

    Observations SBSearchDatabaseSqlite3::find_observations(const string &source, const double mjd_start, double mjd_stop, const int64 limit, const int64 offset) const
    {
        error_if_closed();

        sqlite3_stmt *statement;

        sqlite3_prepare_v2(db, "SELECT observation_id, source, observatory, product_id, mjd_start, mjd_stop, fov, terms FROM observations WHERE source = ? AND mjd_start >= ? AND mjd_stop <= ? LIMIT ? OFFSET ?;", -1, &statement, NULL);
        sqlite3_bind_text(statement, 1, source.c_str(), -1, SQLITE_TRANSIENT);
        sqlite3_bind_double(statement, 2, mjd_start);
        sqlite3_bind_double(statement, 3, mjd_stop);
        sqlite3_bind_int64(statement, 4, limit);
        sqlite3_bind_int64(statement, 5, offset);

        int rc = sqlite3_step(statement);
        check_rc(rc);

        Observations observations;
        observations.reserve(limit);
        while (rc == SQLITE_ROW)
        {
            observations.push_back({string((char *)sqlite3_column_text(statement, 1)),
                                    string((char *)sqlite3_column_text(statement, 2)),
                                    string((char *)sqlite3_column_text(statement, 3)),
                                    sqlite3_column_double(statement, 4),
                                    sqlite3_column_double(statement, 5),
                                    string((char *)sqlite3_column_text(statement, 6)),
                                    string((char *)sqlite3_column_text(statement, 7)),
                                    sqlite3_column_int64(statement, 0)});
            rc = sqlite3_step(statement);
            check_rc(rc);
        }

        sqlite3_finalize(statement);
        return observations;
    }

    Observations SBSearchDatabaseSqlite3::find_observations(vector<string> query_terms, const Options &options) const
    {
        // query_terms may be spatial-temporal, just spatial, or just temporal.
        error_if_closed();

        int rc;
        int count = 0;
        string term_string;
        std::set<int64> approximate_matches;
        sqlite3_stmt *stmt;

        // Query database with terms, but not too many at once
        string statement;
        statement.reserve(MAXIMUM_QUERY_TERMS * 15);
        for (size_t i = 0; i < query_terms.size(); i += MAXIMUM_QUERY_TERMS)
        {
            int chunk = std::min(query_terms.size() - i, MAXIMUM_QUERY_TERMS);
            statement.clear();
            statement.append("SELECT rowid FROM observations_terms_index WHERE terms MATCH '");

            // append quoted terms
            statement += '"' + query_terms[i] + '"';
            // join with OR
            for (size_t j = 1; j < chunk; j++)
            {
                statement.append(" OR \"");
                statement.append(query_terms[i + j]);
                statement.append("\"");
            }

            statement.append("'");

            sqlite3_prepare_v2(db, statement.c_str(), -1, &stmt, NULL);

            sqlite3_bind_double(stmt, options.mjd_start, options.mjd_stop);
            if (!options.source.empty())
                sqlite3_bind_text(stmt, 3, options.source.c_str(), -1, SQLITE_TRANSIENT);

            rc = sqlite3_step(stmt);
            check_rc(rc);
            while (rc == SQLITE_ROW)
            {
                approximate_matches.insert(sqlite3_column_int64(stmt, 0));
                rc = sqlite3_step(stmt);
                check_rc(rc);
            }
            sqlite3_finalize(stmt);
            count += chunk;
        }

        Logger::debug() << "Searched " << count << " of " << query_terms.size() << " query terms."
                        << endl;

        Observations observations = get_observations(approximate_matches.begin(), approximate_matches.end());
        observations.erase(std::remove_if(observations.begin(), observations.end(),
                                          [mjd_start = options.mjd_start, mjd_stop = options.mjd_stop](const Observation &obs)
                                          { return ((obs.mjd_start() < mjd_start) | (obs.mjd_stop() > mjd_stop)); }),
                           observations.end());
        if (!options.source.empty())
            observations.erase(std::remove_if(observations.begin(), observations.end(),
                                              [source = options.source](const Observation &obs)
                                              { return obs.source() != source; }),
                               observations.end());
        return observations;
    }

    void SBSearchDatabaseSqlite3::add_found(const Found &found) const
    {
        error_if_closed();

        char now[32];
        std::time_t time_now = std::time(nullptr);
        std::strftime(now, 32, "%F %T", std::gmtime(&time_now));

        int rc;
        sqlite3_stmt *stmt;
        sqlite3_prepare_v2(db, "INSERT INTO found VALUES (NULL, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);", -1, &stmt, NULL);
        sqlite3_bind_int64(stmt, 1, found.observation.observation_id());
        sqlite3_bind_int64(stmt, 2, found.ephemeris.target().moving_target_id());
        sqlite3_bind_double(stmt, 3, found.ephemeris.data(0).mjd);
        sqlite3_bind_double(stmt, 4, found.ephemeris.data(0).tmtp);
        sqlite3_bind_double(stmt, 5, found.ephemeris.data(0).ra);
        sqlite3_bind_double(stmt, 6, found.ephemeris.data(0).dec);
        sqlite3_bind_double(stmt, 7, found.ephemeris.data(0).unc_a);
        sqlite3_bind_double(stmt, 8, found.ephemeris.data(0).unc_b);
        sqlite3_bind_double(stmt, 9, found.ephemeris.data(0).unc_theta);
        sqlite3_bind_double(stmt, 10, found.ephemeris.data(0).rh);
        sqlite3_bind_double(stmt, 11, found.ephemeris.data(0).delta);
        sqlite3_bind_double(stmt, 12, found.ephemeris.data(0).phase);
        sqlite3_bind_double(stmt, 13, found.ephemeris.data(0).selong);
        sqlite3_bind_double(stmt, 14, found.ephemeris.data(0).true_anomaly);
        sqlite3_bind_double(stmt, 15, found.ephemeris.data(0).sangle);
        sqlite3_bind_double(stmt, 16, found.ephemeris.data(0).vangle);
        sqlite3_bind_double(stmt, 17, found.ephemeris.data(0).vmag);
        sqlite3_bind_text(stmt, 18, now, -1, SQLITE_STATIC);

        rc = sqlite3_step(stmt);
        check_rc(rc);
        sqlite3_finalize(stmt);
    }

    Founds SBSearchDatabaseSqlite3::get_found(const Observation &observation) const
    {
        int rc;
        sqlite3_stmt *stmt;

        sqlite3_prepare_v2(db, R"(
SELECT
    moving_target_id, mjd, tmtp,
    ra, dec, unc_a, unc_b, unc_theta,
    rh, delta, phase, selong, true_anomaly,
    sangle, vangle, vmag
FROM found
WHERE observation_id=?;
)",
                           -1, &stmt, NULL);

        sqlite3_bind_int64(stmt, 1, observation.observation_id());
        rc = sqlite3_step(stmt);
        check_rc(rc);

        Founds founds;
        while (rc == SQLITE_ROW)
        {
            MovingTarget target = get_moving_target(sqlite3_column_int(stmt, 0));

            Ephemeris::Datum d;
            d.mjd = sqlite3_column_double(stmt, 1);
            d.tmtp = sqlite3_column_double(stmt, 2);
            d.ra = sqlite3_column_double(stmt, 3);
            d.dec = sqlite3_column_double(stmt, 4);
            d.unc_a = sqlite3_column_double(stmt, 5);
            d.unc_b = sqlite3_column_double(stmt, 6);
            d.unc_theta = sqlite3_column_double(stmt, 7);
            d.rh = sqlite3_column_double(stmt, 8);
            d.delta = sqlite3_column_double(stmt, 9);
            d.phase = sqlite3_column_double(stmt, 10);
            d.selong = sqlite3_column_double(stmt, 11);
            d.true_anomaly = sqlite3_column_double(stmt, 12);
            d.sangle = sqlite3_column_double(stmt, 13);
            d.vangle = sqlite3_column_double(stmt, 14);
            d.vmag = sqlite3_column_double(stmt, 15);

            Ephemeris ephemeris(target, {d});
            founds.append(Found(observation, ephemeris));

            rc = sqlite3_step(stmt);
            check_rc(rc);
        }

        sqlite3_finalize(stmt);

        return founds;
    }

    Founds SBSearchDatabaseSqlite3::get_found(const MovingTarget &target) const
    {
        int rc;
        sqlite3_stmt *stmt;

        sqlite3_prepare_v2(db, R"(
SELECT
    observation_id, mjd, tmtp,
    ra, dec, unc_a, unc_b, unc_theta,
    rh, delta, phase, selong, true_anomaly,
    sangle, vangle, vmag
FROM found
WHERE moving_target_id=?;
)",
                           -1, &stmt, NULL);

        sqlite3_bind_int64(stmt, 1, target.moving_target_id());
        rc = sqlite3_step(stmt);
        check_rc(rc);

        Founds founds;
        while (rc == SQLITE_ROW)
        {
            Observation observation = get_observation(sqlite3_column_int64(stmt, 0));

            Ephemeris::Datum d;
            d.mjd = sqlite3_column_double(stmt, 1);
            d.tmtp = sqlite3_column_double(stmt, 2);
            d.ra = sqlite3_column_double(stmt, 3);
            d.dec = sqlite3_column_double(stmt, 4);
            d.unc_a = sqlite3_column_double(stmt, 5);
            d.unc_b = sqlite3_column_double(stmt, 6);
            d.unc_theta = sqlite3_column_double(stmt, 7);
            d.rh = sqlite3_column_double(stmt, 8);
            d.delta = sqlite3_column_double(stmt, 9);
            d.phase = sqlite3_column_double(stmt, 10);
            d.selong = sqlite3_column_double(stmt, 11);
            d.true_anomaly = sqlite3_column_double(stmt, 12);
            d.sangle = sqlite3_column_double(stmt, 13);
            d.vangle = sqlite3_column_double(stmt, 14);
            d.vmag = sqlite3_column_double(stmt, 15);

            Ephemeris ephemeris(target, {d});
            founds.append(Found(observation, ephemeris));

            rc = sqlite3_step(stmt);
            check_rc(rc);
        }

        sqlite3_finalize(stmt);

        return founds;
    }

    void SBSearchDatabaseSqlite3::remove_found(const Found &found) const
    {
        // found rows are unique by observation_id and moving_target_id
        int rc;
        sqlite3_stmt *stmt;
        sqlite3_prepare_v2(db, "DELETE FROM found WHERE moving_target_id=? and observation_id=?;", -1, &stmt, NULL);
        sqlite3_bind_int64(stmt, 1, found.ephemeris.target().moving_target_id());
        sqlite3_bind_int64(stmt, 2, found.observation.observation_id());
        rc = sqlite3_step(stmt);
        check_rc(rc);
        sqlite3_finalize(stmt);
    }

    void SBSearchDatabaseSqlite3::check_rc(const int rc) const
    {
        if ((rc != SQLITE_OK) & (rc != SQLITE_ROW) & (rc != SQLITE_DONE))
        {
            Logger::error() << "sqlite3 error (" << rc << ") " << sqlite3_errmsg(db) << endl;
            throw std::runtime_error("sqlite3 error");
        }
    }

    void SBSearchDatabaseSqlite3::check_sql(char *error_message) const
    {
        error_if_closed();

        if (error_message != NULL)
        {
            Logger::error() << error_message << endl;
            sqlite3_free(error_message);
            throw std::runtime_error("\nSQL error\n");
        }
    }

    void SBSearchDatabaseSqlite3::error_if_closed() const
    {
        if (db == NULL)
            throw std::runtime_error("Database is not open.");
    }
}