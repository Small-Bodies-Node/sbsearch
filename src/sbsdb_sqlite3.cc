#include "config.h"

#include <algorithm>
#include <ctime>
#include <iostream>
#include <stdexcept>
#include <set>
#include <string>
#include <vector>
#include <sqlite3.h>

#include "ephemeris.h"
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
using std::to_string;
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

        setup_tables();
        observatories_ = get_observatories();
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
        // Consider a tokenizer that includes hyphens.
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
  object_id INTEGER NOT NULL,
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
)");

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
        vector<string> values = {to_string(options.max_spatial_cells()),
                                 to_string(options.max_spatial_level()),
                                 to_string(options.min_spatial_level()),
                                 to_string(options.temporal_resolution())};
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
            throw MovingTargetError("object_id " + to_string(target.object_id()) + " already exists");
        sqlite3_finalize(stmt);

        Logger::info() << "Add moving target " << target.designation() << endl;
        try
        {
            execute_sql("BEGIN TRANSACTION;");
            add_moving_target_name(object_id, target.designation(), true);
            for (const string &name : target.alternate_names())
                add_moving_target_name(target.object_id(), name, false);
            execute_sql("END TRANSACTION;");
        }
        catch (const MovingTargetError &err)
        {
            execute_sql("ROLLBACK TRANSACTION;");
            throw err;
        }
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
        if (rc == SQLITE_CONSTRAINT)
            throw MovingTargetError("Target name already exists in the database " + name);
        check_rc(rc);
        sqlite3_finalize(stmt);
    };

    void SBSearchDatabaseSqlite3::remove_moving_target(const MovingTarget &target)
    {
        error_if_closed();
        int rc;
        sqlite3_stmt *stmt;

        try
        {
            execute_sql("BEGIN TRANSACTION;");
            sqlite3_prepare_v2(db, "SELECT COUNT() FROM moving_targets WHERE object_id=?", -1, &stmt, NULL);
            sqlite3_bind_int(stmt, 1, target.object_id());
            rc = sqlite3_step(stmt);
            check_rc(rc);
            int count = sqlite3_column_int(stmt, 0);
            if (count == 0)
                throw MovingTargetError("object_id " + to_string(target.object_id()) + " not found");
            sqlite3_finalize(stmt);

            Logger::info() << "Remove moving target " << target.designation() << endl;
            sqlite3_prepare_v2(db, "DELETE FROM moving_targets WHERE object_id=?", -1, &stmt, NULL);
            sqlite3_bind_int(stmt, 1, target.object_id());
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
    }

    void SBSearchDatabaseSqlite3::update_moving_target(const MovingTarget &target)
    {
        Logger::info() << "Update moving target " << target.designation() << endl;
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
            throw MovingTargetError("object_id " + to_string(target.object_id()) + " not found");

        // otherwise, loop through the names and add them to our object
        while (rc == SQLITE_ROW)
        {
            bool primary = (bool)sqlite3_column_int(stmt, 1);
            string name(reinterpret_cast<const char *>(sqlite3_column_text(stmt, 0)));

            if (primary)
                target.designation(name);
            else
                target.add_name(name);

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

    void SBSearchDatabaseSqlite3::add_observatory(const string &name, const Observatory &observatory)
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

        observatories_[name] = observatory;
    }

    const Observatory SBSearchDatabaseSqlite3::get_observatory(const string &name)
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
            throw ObservatoryError(name + "not found");

        Observatory observatory{sqlite3_column_double(stmt, 0),
                                sqlite3_column_double(stmt, 1),
                                sqlite3_column_double(stmt, 2)};
        sqlite3_finalize(stmt);
        return observatory;
    }

    const Observatories SBSearchDatabaseSqlite3::get_observatories()
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

        observatories_ = observatories;
        return observatories;
    }

    void SBSearchDatabaseSqlite3::remove_observatory(const string &name)
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

        observatories_.erase(name);
    }

    void SBSearchDatabaseSqlite3::add_ephemeris(Ephemeris &eph)
    {
        error_if_closed();

        // validate ephemeris target
        if (eph.target().object_id() == UNDEF_OBJECT_ID)
        {
            MovingTarget target = eph.target();
            add_moving_target(target);
            eph.target(target); // update ephemeris object
        }
        else
        {
            MovingTarget target;
            try
            {
                target = get_moving_target(eph.target().object_id());
                if (eph.target() != target)
                    throw MovingTargetError("Ephemeris target does not match database copy.");
            }
            catch (const MovingTargetError &error)
            {
                // object ID was not in the database, so we can add this as a new target
                target = eph.target();
                add_moving_target(target);
                eph.target(target); // update ephemeris object
            }
        }

        Logger::info() << "Adding " << to_string(eph.num_vertices()) << " ephemeris epochs for target " << eph.target().designation() << " (object_id=" << eph.target().object_id() << ")." << endl;

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
  object_id, mjd, tmtp,
  ra, dec, unc_a, unc_b, unc_theta,
  rh, delta, phase, selong, true_anomaly,
  sangle, vangle, vmag, retrieved
) VALUES (
  ?, ?, ?, ?, ?, ?, ?, ?,
  ?, ?, ?, ?, ?, ?, ?, ?, ?
);)",
                               -1, &stmt, NULL);
            sqlite3_bind_int(stmt, 1, eph.target().object_id());
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

    Ephemeris SBSearchDatabaseSqlite3::get_ephemeris(const MovingTarget target, double mjd_start, double mjd_stop)
    {
        int rc;
        sqlite3_stmt *stmt;
        Ephemeris::Data data;

        sqlite3_prepare_v2(db, "SELECT COUNT() FROM ephemerides WHERE object_id=? AND mjd >= ? and mjd <= ?;", -1, &stmt, NULL);
        sqlite3_bind_int(stmt, 1, target.object_id());
        sqlite3_bind_double(stmt, 2, mjd_start);
        sqlite3_bind_double(stmt, 3, mjd_stop);
        rc = sqlite3_step(stmt);
        check_rc(rc);
        int count = sqlite3_column_int(stmt, 0);
        sqlite3_finalize(stmt);
        Logger::debug() << "Reading " << count << " ephemeris epochs from database for " << target.designation() << " (object_id=" << target.object_id() << ")." << endl;

        data.reserve(count);
        sqlite3_prepare_v2(db, R"(
SELECT
    mjd, tmtp, ra, dec, unc_a, unc_b, unc_theta,
    rh, delta, phase, selong, true_anomaly,
    sangle, vangle, vmag
FROM ephemerides
WHERE object_id=? AND mjd >= ? and mjd <= ?;)",
                           -1, &stmt, NULL);
        sqlite3_bind_int(stmt, 1, target.object_id());
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

    int SBSearchDatabaseSqlite3::remove_ephemeris(const MovingTarget target, double mjd_start, double mjd_stop)
    {
        int rc;
        sqlite3_stmt *stmt;

        sqlite3_prepare_v2(db, "SELECT COUNT() FROM ephemerides WHERE object_id=? AND mjd >= ? and mjd <= ?;", -1, &stmt, NULL);
        sqlite3_bind_int(stmt, 1, target.object_id());
        sqlite3_bind_double(stmt, 2, mjd_start);
        sqlite3_bind_double(stmt, 3, mjd_stop);
        rc = sqlite3_step(stmt);
        check_rc(rc);
        int count = sqlite3_column_int(stmt, 0);
        sqlite3_finalize(stmt);
        Logger::info() << "Removing " << count << " ephemeris epochs from database for " << target.designation() << " (object_id=" << target.object_id() << ")." << endl;

        sqlite3_prepare_v2(db, "DELETE FROM ephemerides WHERE object_id=? AND mjd >= ? and mjd <= ?;",
                           -1, &stmt, NULL);
        sqlite3_bind_int(stmt, 1, target.object_id());
        sqlite3_bind_double(stmt, 2, mjd_start);
        sqlite3_bind_double(stmt, 3, mjd_stop);
        rc = sqlite3_step(stmt);
        check_rc(rc);
        sqlite3_finalize(stmt);

        return count;
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
            rc = sqlite3_prepare_v2(db, "INSERT INTO observations VALUES (NULL, ?, ?, ?, ?, ?, ?, ?) RETURNING observation_id;", -1, &statement, NULL);
            check_rc(rc);
        }
        else
        {
            // update existing observation
            rc = sqlite3_prepare_v2(db, "UPDATE observations SET source=?, observatory=?, product_id=?, mjd_start=?, mjd_stop=?, fov=?, terms=? WHERE observation_id=?;", -1, &statement, NULL);
            rc += sqlite3_bind_int64(statement, 8, observation.observation_id());
        }

        rc += sqlite3_bind_text(statement, 1, observation.source().c_str(), -1, SQLITE_TRANSIENT);
        rc += sqlite3_bind_text(statement, 2, observation.observatory().c_str(), -1, SQLITE_TRANSIENT);
        rc += sqlite3_bind_text(statement, 3, observation.product_id().c_str(), -1, SQLITE_TRANSIENT);
        rc += sqlite3_bind_double(statement, 4, observation.mjd_start());
        rc += sqlite3_bind_double(statement, 5, observation.mjd_stop());
        rc += sqlite3_bind_text(statement, 6, observation.fov().c_str(), -1, SQLITE_TRANSIENT);
        rc += sqlite3_bind_text(statement, 7, observation.terms().c_str(), -1, SQLITE_TRANSIENT);

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

    vector<Observation> SBSearchDatabaseSqlite3::find_observations(const double mjd_start, double mjd_stop)
    {
        error_if_closed();

        sqlite3_stmt *statement;

        sqlite3_prepare_v2(db, "SELECT observation_id, source, observatory, product_id, mjd_start, mjd_stop, fov, terms FROM observations WHERE mjd_start > ? AND mjd_stop < ?;", -1, &statement, NULL);
        sqlite3_bind_double(statement, 1, mjd_start);
        sqlite3_bind_double(statement, 2, mjd_stop);

        int rc = sqlite3_step(statement);
        check_rc(rc);

        vector<Observation> observations;
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

    vector<Observation> SBSearchDatabaseSqlite3::find_observations(const string &source, const double mjd_start, double mjd_stop)
    {
        error_if_closed();

        sqlite3_stmt *statement;

        sqlite3_prepare_v2(db, "SELECT observation_id, source, observatory, product_id, mjd_start, mjd_stop, fov, terms FROM observations WHERE source = ? AND mjd_start > ? AND mjd_stop < ?;", -1, &statement, NULL);
        sqlite3_bind_text(statement, 1, source.c_str(), -1, SQLITE_TRANSIENT);
        sqlite3_bind_double(statement, 2, mjd_start);
        sqlite3_bind_double(statement, 3, mjd_stop);

        int rc = sqlite3_step(statement);
        check_rc(rc);

        vector<Observation> observations;
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

    vector<Observation> SBSearchDatabaseSqlite3::find_observations(vector<string> query_terms, const Options &options)
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

        vector<Observation> observations = get_observations(approximate_matches.begin(), approximate_matches.end());
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