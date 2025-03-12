#include "config.h"

#include <algorithm>
#include <ctime>
#include <iostream>
#include <stdexcept>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <pqxx/pqxx>

#include "ephemeris.h"
#include "found.h"
#include "exceptions.h"
#include "logging.h"
#include "moving_target.h"
#include "observation.h"
#include "observatory.h"
#include "sbsdb.h"
#include "sbsdb_postgresql.h"

using std::endl;
using std::set;
using std::string;
using std::vector;

namespace sbsearch
{
    void SBSearchDatabasePostgreSQL::close()
    {
        if (connection_.is_open())
        {
            connection_.close();
            Logger::info() << "Closed database." << endl;
        }
    };

    void SBSearchDatabasePostgreSQL::execute_sql(const char *statement)
    {
        error_if_closed();

        try
        {
            pqxx::nontransaction work(connection_);
            work.exec(statement);
        }
        catch (std::exception const &e)
        {
            std::cerr << e.what() << std::endl;
            return;
        }
    };

    void SBSearchDatabasePostgreSQL::setup_tables()
    {
        error_if_closed();

        {
            pqxx::nontransaction work(connection_);

            work.exec(R"(
CREATE TABLE IF NOT EXISTS observations (
  observation_id SERIAL PRIMARY KEY,
  source TEXT NOT NULL,
  observatory VARCHAR(64) NOT NULL,
  product_id VARCHAR(128) NOT NULL,
  mjd_start DOUBLE PRECISION NOT NULL,
  mjd_stop DOUBLE PRECISION NOT NULL,
  fov VARCHAR(128) NOT NULL,
  terms TEXT[] NOT NULL
);
)");
        }

        create_observations_indices();

        {
            pqxx::nontransaction work(connection_);
            work.exec(R"(
CREATE TABLE IF NOT EXISTS moving_targets (
  moving_targets_row_id SERIAL PRIMARY KEY,
  moving_target_id INTEGER NOT NULL,
  name VARCHAR(64) NOT NULL,
  small_body BOOLEAN NOT NULL,
  primary_id BOOLEAN NOT NULL
);
CREATE UNIQUE INDEX IF NOT EXISTS idx_moving_target_primary_id ON moving_targets(moving_target_id, primary_id) WHERE primary_id=TRUE;
CREATE UNIQUE INDEX IF NOT EXISTS idx_moving_target_name_small_body ON moving_targets(name, small_body);
CREATE INDEX IF NOT EXISTS idx_moving_target_moving_target_id ON moving_targets(moving_target_id);

CREATE TABLE IF NOT EXISTS observatories (
  name VARCHAR(64) UNIQUE NOT NULL,
  longitude DOUBLE PRECISION NOT NULL,
  rho_cos_phi DOUBLE PRECISION NOT NULL,
  rho_sin_phi DOUBLE PRECISION NOT NULL
);

CREATE TABLE IF NOT EXISTS ephemerides (
  ephemeris_id SERIAL PRIMARY KEY,
  moving_target_id INTEGER NOT NULL,
  mjd DOUBLE PRECISION NOT NULL,
  tmtp DOUBLE PRECISION NOT NULL,
  ra DOUBLE PRECISION NOT NULL,
  dec DOUBLE PRECISION NOT NULL,
  unc_a DOUBLE PRECISION NOT NULL,
  unc_b DOUBLE PRECISION NOT NULL,
  unc_theta DOUBLE PRECISION NOT NULL,
  rh DOUBLE PRECISION NOT NULL,
  delta DOUBLE PRECISION NOT NULL,
  phase DOUBLE PRECISION NOT NULL,
  selong DOUBLE PRECISION NOT NULL,
  true_anomaly DOUBLE PRECISION NOT NULL,
  sangle DOUBLE PRECISION NOT NULL,
  vangle DOUBLE PRECISION NOT NULL,
  vmag DOUBLE PRECISION NOT NULL,
  retrieved VARCHAR(64) NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_ephemerides_moving_target_id ON ephemerides(moving_target_id);

CREATE TABLE IF NOT EXISTS found (
  found_id SERIAL PRIMARY KEY,
  observation_id INTEGER NOT NULL,
  moving_target_id INTEGER NOT NULL,
  mjd DOUBLE PRECISION NOT NULL,
  tmtp DOUBLE PRECISION NOT NULL,
  ra DOUBLE PRECISION NOT NULL,
  dec DOUBLE PRECISION NOT NULL,
  unc_a DOUBLE PRECISION NOT NULL,
  unc_b DOUBLE PRECISION NOT NULL,
  unc_theta DOUBLE PRECISION NOT NULL,
  rh DOUBLE PRECISION NOT NULL,
  delta DOUBLE PRECISION NOT NULL,
  phase DOUBLE PRECISION NOT NULL,
  selong DOUBLE PRECISION NOT NULL,
  true_anomaly DOUBLE PRECISION NOT NULL,
  sangle DOUBLE PRECISION NOT NULL,
  vangle DOUBLE PRECISION NOT NULL,
  vmag DOUBLE PRECISION,
  saved VARCHAR(64) NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_found_observation_id ON found(observation_id);
CREATE UNIQUE INDEX IF NOT EXISTS idx_found_moving_target_id ON found(moving_target_id, observation_id);

CREATE TABLE IF NOT EXISTS configuration (
  parameter VARCHAR(64) NOT NULL UNIQUE,
  value TEXT NOT NULL
);
INSERT INTO configuration VALUES ('max_spatial_index_cells', '8') ON CONFLICT DO NOTHING;
INSERT INTO configuration VALUES ('max_spatial_level', '12') ON CONFLICT DO NOTHING;
INSERT INTO configuration VALUES ('min_spatial_level', '4') ON CONFLICT DO NOTHING;
INSERT INTO configuration VALUES ('temporal_resolution', '1') ON CONFLICT DO NOTHING;
INSERT INTO configuration VALUES ('database version', ')" SBSEARCH_DATABASE_VERSION R"(') ON CONFLICT DO NOTHING;

ANALYZE;
)");
        }
        Logger::debug() << "Database tables are set." << endl;
    };

    void SBSearchDatabasePostgreSQL::drop_observations_indices()
    {
        Logger::info() << "Dropping observations indices." << std::endl;
        execute_sql("DROP INDEX idx_observations_terms;");
        execute_sql("DROP INDEX idx_observations_mjd_start;");
        execute_sql("DROP INDEX idx_observations_mjd_stop;");
        execute_sql("DROP INDEX idx_observations_source_mjd_start;");
        execute_sql("DROP INDEX idx_observations_source_mjd_stop;");
        execute_sql("DROP INDEX idx_observations_observatory;");
        execute_sql("DROP INDEX idx_observations_product_id;");
        Logger::info() << "Observations indices dropped." << std::endl;
    };

    void SBSearchDatabasePostgreSQL::create_observations_indices()
    {
        Logger::info() << "Creating observations indices." << std::endl;

        execute_sql(R"(
        CREATE INDEX IF NOT EXISTS idx_observations_terms
        ON observations
        USING GIN (terms);

        CREATE INDEX IF NOT EXISTS idx_observations_mjd_start
        ON observations(mjd_start);

        CREATE INDEX IF NOT EXISTS idx_observations_mjd_stop
        ON observations(mjd_stop);

        CREATE INDEX IF NOT EXISTS idx_observations_source_mjd_start
        ON observations(source, mjd_start);

        CREATE INDEX IF NOT EXISTS idx_observations_source_mjd_stop
        ON observations(source, mjd_stop);

        CREATE INDEX IF NOT EXISTS idx_observations_observatory
        ON observations(observatory);

        CREATE UNIQUE INDEX IF NOT EXISTS idx_observations_product_id
        ON observations(product_id);

        ANALYZE observations;
)");

        Logger::info()
            << "Created observations indices." << std::endl;
    };

    // get single value results from a SQL statement
    optional<double> SBSearchDatabasePostgreSQL::get_double(const char *statement)
    {
        error_if_closed();

        pqxx::nontransaction work(connection_);
        pqxx::row row = work.exec1(statement);
        if (row[0].is_null())
            return {};

        return {row[0].as<double>()};
    };

    optional<int> SBSearchDatabasePostgreSQL::get_int(const char *statement)
    {
        error_if_closed();
        pqxx::nontransaction work(connection_);
        pqxx::row row = work.exec1(statement);
        if (row[0].is_null())
            return {};

        return {row[0].as<int>()};
    };

    optional<int64_t> SBSearchDatabasePostgreSQL::get_int64(const char *statement)
    {
        error_if_closed();
        pqxx::nontransaction work(connection_);
        pqxx::row row = work.exec1(statement);
        if (row[0].is_null())
            return {};

        return {row[0].as<int64_t>()};
    };

    optional<string> SBSearchDatabasePostgreSQL::get_string(const char *statement)
    {
        error_if_closed();
        pqxx::nontransaction work(connection_);
        pqxx::row row = work.exec1(statement);
        if (row[0].is_null())
            return {};

        return {row[0].as<string>()};
    };

    void SBSearchDatabasePostgreSQL::indexer_options(Indexer::Options options)
    {
        error_if_closed();

        vector<string> parameters = {"max_spatial_index_cells",
                                     "max_spatial_level",
                                     "min_spatial_level",
                                     "temporal_resolution"};
        vector<string> values = {std::to_string(options.max_spatial_index_cells()),
                                 std::to_string(options.max_spatial_level()),
                                 std::to_string(options.min_spatial_level()),
                                 std::to_string(options.temporal_resolution())};

        pqxx::work work(connection_);
        for (int i = 0; i < parameters.size(); i++)
            work.exec_params("UPDATE configuration SET value = $1 WHERE parameter = $2", values[i], parameters[i]);

        work.commit();
    }

    std::pair<optional<double>, optional<double>> SBSearchDatabasePostgreSQL::observation_date_range(const string &source)
    {
        error_if_closed();

        optional<double> mjd_start;
        optional<double> mjd_stop;

        if (source.empty() | (source == ""))
        {
            mjd_start = get_double("SELECT MIN(mjd_start) FROM observations");
            mjd_stop = get_double("SELECT MAX(mjd_stop) FROM observations");
        }
        else
        {
            pqxx::nontransaction work(connection_);
            pqxx::row row = work.exec_params1(
                "SELECT MIN(mjd_start), MAX(mjd_stop) FROM observations WHERE source=$1",
                source);

            if (!row[0].is_null())
                mjd_start = row[0].as<double>();

            if (!row[1].is_null())
                mjd_stop = row[1].as<double>();
        }

        return {mjd_start, mjd_stop};
    };

    void SBSearchDatabasePostgreSQL::add_moving_target(MovingTarget &target)
    {
        error_if_closed();

        if (!target.moving_target_id())
        {
            // new objects use db largest moving_target_id + 1, or 1 if there are no objects
            // target.moving_target_id(get_int("SELECT IFNULL(MAX(moving_target_id), 0) + 1 FROM moving_targets"));
            target.moving_target_id(
                get_int("SELECT COALESCE(MAX(moving_target_id), 0) + 1 "
                        "FROM moving_targets"));
        }

        pqxx::nontransaction work(connection_);
        pqxx::row row = work.exec_params1(
            "SELECT COUNT(*) FROM moving_targets WHERE moving_target_id=$1",
            target.moving_target_id().value());
        int count = row[0].as<int>();
        if (count != 0)
            throw MovingTargetError("moving target id " +
                                    std::to_string(target.moving_target_id().value()) +
                                    " already exists");

        Logger::info() << "Add moving target " << target << endl;
        try
        {
            work.exec("BEGIN TRANSACTION;");
            add_moving_target_name(work, target.moving_target_id().value(), target.designation(), target.small_body(), true);
            for (const string &name : target.alternate_names())
                add_moving_target_name(work, target.moving_target_id().value(), name, target.small_body(), false);
            work.exec("END TRANSACTION;");
        }
        catch (const MovingTargetError &err)
        {
            work.exec("ROLLBACK TRANSACTION;");
            throw err;
        }
        Logger::info() << target << " added to database." << std::endl;
    };

    void SBSearchDatabasePostgreSQL::add_moving_target_name(pqxx::transaction_base &work,
                                                            const int64_t moving_target_id,
                                                            const string &name,
                                                            const bool small_body,
                                                            const bool primary_id)
    {
        Logger::debug() << "Add moving target " << name
                        << " (ID=" << moving_target_id
                        << "; small body=" << (small_body ? "true" : "false")
                        << "; primary=" << (primary_id ? "true" : "false")
                        << ")" << std::endl;

        error_if_closed();

        try
        {
            work.exec_params(
                "INSERT INTO moving_targets (moving_target_id, name, small_body, primary_id) "
                "VALUES ($1, $2, $3, $4)",
                moving_target_id,
                name,
                small_body,
                primary_id);
        }
        catch (const pqxx::unique_violation &e)
        {
            throw MovingTargetError("target uniqueness violated: " + string(e.what()));
        }
    };

    void SBSearchDatabasePostgreSQL::remove_moving_target(const MovingTarget &target)
    {
        error_if_closed();

        pqxx::nontransaction work(connection_);

        try
        {
            work.exec("BEGIN TRANSACTION");
            pqxx::row row = work.exec_params1(
                "SELECT COUNT(*) FROM moving_targets "
                "WHERE moving_target_id=$1 AND small_body=$2",
                target.moving_target_id().value(),
                target.small_body());

            if (row[0].as<int>() == 0)
                throw MovingTargetError("moving target id " +
                                        std::to_string(target.moving_target_id().value()) +
                                        " with small body flag " +
                                        (target.small_body() ? "true" : "false") +
                                        " not found");

            Logger::info() << "Remove moving target " << target << endl;
            work.exec_params0(
                "DELETE FROM moving_targets WHERE moving_target_id=$1 AND small_body=$2",
                target.moving_target_id().value(),
                target.small_body());
            work.exec("END TRANSACTION");
        }
        catch (const MovingTargetError &err)
        {
            work.exec("ROLLBACK TRANSACTION");
            throw err;
        }
        Logger::info() << target << " removed from database." << std::endl;
    };

    MovingTarget SBSearchDatabasePostgreSQL::get_moving_target(const int64_t moving_target_id)
    {
        error_if_closed();
        MovingTarget target;
        target.moving_target_id(moving_target_id);

        pqxx::nontransaction work(connection_);

        pqxx::result result = work.exec_params(
            "SELECT name, small_body, primary_id FROM moving_targets WHERE moving_target_id=$1",
            moving_target_id);

        if (result.size() == 0)
            throw MovingTargetError("moving target id " +
                                    std::to_string(target.moving_target_id().value()) +
                                    " not found");

        target.small_body(result[0][1].as<bool>());

        // loop through the names and add them to our object
        for (auto const &row : result)
        {
            bool primary = row[2].as<bool>();
            string name = row[0].as<string>();

            if (primary)
                target.designation(name);
            else
                target.add_name(name);
        }

        return target;
    };

    MovingTarget SBSearchDatabasePostgreSQL::get_moving_target(const string &name, const bool small_body)
    {
        error_if_closed();
        int64_t moving_target_id;

        try
        {
            pqxx::nontransaction work(connection_);
            pqxx::row row = work.exec_params1(
                "SELECT moving_target_id FROM moving_targets "
                "WHERE name=$1 AND small_body=$2 LIMIT 1",
                name,
                small_body);
            moving_target_id = row[0].as<int64_t>();
        }
        catch (pqxx::unexpected_rows)
        {
            // this name-small body / name combo is not in the database, return a new object
            return MovingTarget(name, small_body);
        }

        // return the target based on moving_target_id
        return get_moving_target(moving_target_id);
    };

    vector<MovingTarget> SBSearchDatabasePostgreSQL::get_all_moving_targets()
    {
        error_if_closed();
        pqxx::nontransaction work(connection_);
        pqxx::result result = work.exec("SELECT DISTINCT(moving_target_id) FROM moving_targets");

        vector<MovingTarget> targets;
        for (auto const &row : result)
            targets.push_back(get_moving_target(row[0].as<int>()));

        return targets;
    };

    void SBSearchDatabasePostgreSQL::add_observatory(const string &name, const Observatory &observatory)
    {
        Logger::info() << "Adding observatory " << name << "." << std::endl;
        error_if_closed();

        // do not add anything if this name is already in the database
        try
        {
            get_observatory(name);
        }
        catch (const ObservatoryError &e)
        {
            pqxx::nontransaction work(connection_);
            work.exec_params(
                "INSERT INTO observatories (name, longitude, rho_cos_phi, rho_sin_phi) "
                "VALUES ($1, $2, $3, $4)",
                name,
                observatory.longitude,
                observatory.rho_cos_phi,
                observatory.rho_sin_phi);
            return;
        }

        throw ObservatoryError(name + " already exists.");
    };

    const Observatory SBSearchDatabasePostgreSQL::get_observatory(const string &name)
    {
        error_if_closed();
        pqxx::nontransaction work(connection_);

        try
        {
            pqxx::row row = work.exec_params1(
                "SELECT longitude, rho_cos_phi, rho_sin_phi "
                "FROM observatories "
                "WHERE name = $1",
                name);
            Observatory observatory{row[0].as<double>(), row[1].as<double>(), row[2].as<double>()};
            return observatory;
        }
        catch (pqxx::unexpected_rows)
        {
            throw ObservatoryError(name + " not found");
        }
    };

    const Observatories SBSearchDatabasePostgreSQL::get_observatories()
    {
        error_if_closed();
        pqxx::nontransaction work(connection_);
        pqxx::result result = work.exec("SELECT name, longitude, rho_cos_phi, rho_sin_phi FROM observatories");

        Observatories observatories;
        for (auto const &row : result)
        {
            const string name = row[0].as<string>();
            Observatory observatory{row[1].as<double>(), row[2].as<double>(), row[3].as<double>()};
            observatories[name] = observatory;
        }

        return observatories;
    };

    void SBSearchDatabasePostgreSQL::remove_observatory(const string &name)
    {

        error_if_closed();
        pqxx::work work(connection_);

        Logger::info() << "Removing observatory with name " << name << endl;
        work.exec_params0("DELETE FROM observatories WHERE name=$1", name);
        work.commit();
    };

    const vector<string> SBSearchDatabasePostgreSQL::get_sources()
    {
        error_if_closed();
        pqxx::nontransaction work(connection_);
        pqxx::result result = work.exec("SELECT DISTINCT(source) FROM observations");

        vector<string> sources;
        for (auto const &row : result)
            sources.push_back(row[0].as<string>());

        return sources;
    };

    void SBSearchDatabasePostgreSQL::add_ephemeris(Ephemeris &eph)
    {
        Logger::info()
            << "Adding " << std::to_string(eph.num_vertices())
            << " ephemeris epochs for target " << eph.target().designation()
            << " (moving_target_id=" << eph.target().moving_target_id().value_or(-1) << ")." << endl;

        error_if_closed();

        if (!eph.target().moving_target_id())
            throw MovingTargetError("Ephemeris target is not in the database.");

        // verify that the moving target ID exists in the database
        MovingTarget target = get_moving_target(
            eph.target().moving_target_id().value()); // throws MovingTargetError if not found
        if (target != eph.target())
            throw MovingTargetError("Ephemeris target does not match database copy.");

        char now[32];
        std::time_t time_now = std::time(nullptr);
        std::strftime(now, 32, "%F %T", std::gmtime(&time_now));

        pqxx::work work(connection_);
        connection_.prepare("", R"(
            INSERT INTO ephemerides (
                moving_target_id, mjd, tmtp,
                ra, dec, unc_a, unc_b, unc_theta,
                rh, delta, phase, selong, true_anomaly,
                sangle, vangle, vmag, retrieved
            ) VALUES (
                $1, $2, $3,
                $4, $5, $6, $7, $8,
                $9, $10, $11, $12, $13,
                $14, $15, $16, $17
            )
        )");
        for (const Ephemeris::Datum row : eph.data())
            work.exec_prepared(
                "",
                eph.target().moving_target_id(), row.mjd, row.tmtp,
                row.ra, row.dec, row.unc_a, row.unc_b, row.unc_theta,
                row.rh, row.delta, row.phase, row.selong, row.true_anomaly,
                row.sangle, row.vangle, row.vmag, now);
        work.commit();
    };

    Ephemeris SBSearchDatabasePostgreSQL::get_ephemeris(const MovingTarget target, double mjd_start, double mjd_stop)
    {
        error_if_closed();

        pqxx::nontransaction work(connection_);
        pqxx::row row = work.exec_params1(
            "SELECT COUNT(*) FROM ephemerides "
            "WHERE moving_target_id=$1 AND mjd >= $2 and mjd <= $3",
            target.moving_target_id().value(),
            mjd_start,
            mjd_stop);
        int count = row[0].as<int>();

        Logger::debug() << "Reading " << count
                        << " ephemeris epochs from database for " << target.designation()
                        << " (moving_target_id=" << target.moving_target_id().value() << ")"
                        << endl;

        Ephemeris::Data data;
        data.reserve(count);
        pqxx::result result = work.exec_params(
            "SELECT"
            "    mjd, tmtp, ra, dec, unc_a, unc_b, unc_theta,"
            "    rh, delta, phase, selong, true_anomaly,"
            "    sangle, vangle, vmag "
            "FROM ephemerides "
            "WHERE moving_target_id = $1 AND mjd >= $2 and mjd <= $3",
            target.moving_target_id().value(),
            mjd_start,
            mjd_stop);

        for (auto const &row : result)
        {
            Ephemeris::Datum d;
            d.mjd = row[0].as<double>();
            d.tmtp = row[1].as<double>();
            d.ra = row[2].as<double>();
            d.dec = row[3].as<double>();
            d.unc_a = row[4].as<double>();
            d.unc_b = row[5].as<double>();
            d.unc_theta = row[6].as<double>();
            d.rh = row[7].as<double>();
            d.delta = row[8].as<double>();
            d.phase = row[9].as<double>();
            d.selong = row[10].as<double>();
            d.true_anomaly = row[11].as<double>();
            d.sangle = row[12].as<double>();
            d.vangle = row[13].as<double>();
            d.vmag = row[14].as<double>();
            data.push_back(d);
        }

        return {target, data};
    };

    int SBSearchDatabasePostgreSQL::remove_ephemeris(const MovingTarget target, double mjd_start, double mjd_stop)
    {
        error_if_closed();
        pqxx::work work(connection_);

        pqxx::row row = work.exec_params1(
            "SELECT COUNT(*) FROM ephemerides "
            "WHERE moving_target_id = $1 AND mjd >= $2 AND mjd <= $3",
            target.moving_target_id().value(),
            mjd_start, mjd_stop);

        int count = row[0].as<int>();

        Logger::info() << "Removing " << count
                       << " ephemeris epochs from database for " << target.designation()
                       << " (moving_target_id=" << target.moving_target_id().value() << ")." << endl;

        try
        {
            work.exec_params(
                "DELETE FROM ephemerides WHERE moving_target_id = $1 AND mjd >= $2 and mjd <= $3",
                target.moving_target_id(),
                mjd_start,
                mjd_stop);
            work.commit();
        }
        catch (std::exception const &err)
        {
            Logger::error() << err.what() << endl;
            Logger::error() << "Error deleting from ephemerides" << endl;
        }

        return count;
    };

    void SBSearchDatabasePostgreSQL::add_new_observation(pqxx::work &work, Observation &observation)
    {
        // insert row and update observation object with observation_id
        pqxx::row row = work.exec_params1(
            "INSERT INTO observations "
            "(source, observatory, product_id, mjd_start, mjd_stop, fov, terms) "
            "VALUES ($1, $2, $3, $4, $5, $6, $7) "
            "RETURNING observation_id",
            observation.source(),
            observation.observatory(),
            observation.product_id(),
            observation.mjd_start(),
            observation.mjd_stop(),
            observation.fov(),
            observation.terms());
        observation.observation_id(row[0].as<int64_t>());
    }

    void SBSearchDatabasePostgreSQL::update_observation(pqxx::work &work, const Observation &observation)
    {
        // update existing observation
        work.exec_params(
            R"(
                INSERT INTO observations
                VALUES ($1, $2, $3, $4, $5, $6, $7, $8)
                ON CONFLICT (observation_id) DO UPDATE SET
                    source=excluded.source,
                    observatory=excluded.observatory,
                    product_id=excluded.product_id,
                    mjd_start=excluded.mjd_start,
                    mjd_stop=excluded.mjd_stop,
                    fov=excluded.fov,
                    terms=excluded.terms
            )",
            observation.observation_id(),
            observation.source(),
            observation.observatory(),
            observation.product_id(),
            observation.mjd_start(),
            observation.mjd_stop(),
            observation.fov(),
            observation.terms());
    }

    void SBSearchDatabasePostgreSQL::add_observations(Observations &observations)
    {
        Logger::info() << "Adding or updating " << observations.size() << " observation"
                       << (observations.size() == 1 ? "" : "s") << "." << endl;

        error_if_closed();
        pqxx::work work(connection_);
        int added = 0, updated = 0;
        for (vector<Observation>::iterator it = observations.begin(); it < observations.end(); it++)
        {
            if (it->terms().size() == 0)
                throw std::runtime_error("Observation is missing index terms.");

            if (it->observation_id() == UNDEFINED_OBSID)
            {
                add_new_observation(work, *it);
                added++;
            }
            else
            {
                update_observation(work, *it);
                updated++;
            }
        }
        work.commit();
        Logger::info() << "Added " << added << " and updated " << updated << " observations." << endl;
    }

    Observation SBSearchDatabasePostgreSQL::get_observation(const int64_t observation_id)
    {
        error_if_closed();
        pqxx::nontransaction work(connection_);
        try
        {
            pqxx::row row = work.exec_params1(
                "SELECT source, observatory, product_id, mjd_start, mjd_stop, fov, terms "
                "FROM observations WHERE observation_id = $1",
                observation_id);

            string source = row[0].as<string>();
            string observatory = row[1].as<string>();
            string product_id = row[2].as<string>();
            double mjd_start = row[3].as<double>();
            double mjd_stop = row[4].as<double>();
            string fov = row[5].as<string>();

            vector<string> terms;
            auto parsed = row[6].as_array();
            std::pair<pqxx::array_parser::juncture, string> next;
            do
            {
                next = parsed.get_next();
                if (next.first == pqxx::array_parser::juncture::string_value)
                    terms.push_back(next.second);
            } while (next.first != pqxx::array_parser::juncture::done);

            return Observation(source, observatory, product_id, mjd_start, mjd_stop, fov, terms, observation_id);
        }
        catch (const std::exception &err)
        {
            Logger::error() << err.what() << endl;
            throw std::runtime_error("No matching observation");
        }
    };

    void SBSearchDatabasePostgreSQL::remove_observations(const double mjd_start, const double mjd_stop)
    {
        error_if_closed();
        pqxx::work work(connection_);
        work.exec_params(
            "DELETE FROM observations WHERE mjd_start >= $1 AND mjd_stop <= $2",
            mjd_start,
            mjd_stop);
        work.commit();
    }

    void SBSearchDatabasePostgreSQL::remove_observations(const string &source, const double mjd_start, const double mjd_stop)
    {
        error_if_closed();
        pqxx::work work(connection_);
        work.exec_params(
            "DELETE FROM observations WHERE source = $1 AND mjd_start >= $2 AND mjd_stop <= $3",
            source,
            mjd_start,
            mjd_stop);
        work.commit();
    };

    int64_t SBSearchDatabasePostgreSQL::count_observations(const double mjd_start, const double mjd_stop)
    {
        error_if_closed();
        pqxx::nontransaction work(connection_);
        return work.exec_params1(
                       "SELECT COUNT(*) FROM observations WHERE mjd_start >= $1 AND mjd_stop <= $2",
                       mjd_start,
                       mjd_stop)[0]
            .as<int64_t>();
    };

    int64_t SBSearchDatabasePostgreSQL::count_observations(const string &source, const double mjd_start, const double mjd_stop)
    {
        if (source == "")
            return count_observations(mjd_start, mjd_stop);

        error_if_closed();
        pqxx::nontransaction work(connection_);
        return work.exec_prepared1(
                       "SELECT COUNT(*) FROM observations WHERE source = $1 AND mjd_start >= $2 AND mjd_stop <= $3",
                       source,
                       mjd_start,
                       mjd_stop)[0]
            .as<int64_t>();
    };

    Observations SBSearchDatabasePostgreSQL::find_observations(const double mjd_start, const double mjd_stop, const int64_t limit, const int64_t offset)
    {
        error_if_closed();
        pqxx::nontransaction work(connection_);

        pqxx::result result = work.exec_params(
            "SELECT observation_id, source, observatory, product_id, mjd_start, mjd_stop, fov, terms "
            "FROM observations "
            "WHERE mjd_start >= $1 AND mjd_stop <= $2 "
            "LIMIT $3 OFFSET $4",
            mjd_start,
            mjd_stop,
            limit,
            offset);

        Observations observations;
        observations.data.reserve(limit);
        for (auto const &row : result)
            observations.append({row[0].as<string>(),
                                 row[1].as<string>(),
                                 row[2].as<string>(),
                                 row[3].as<double>(),
                                 row[4].as<double>(),
                                 row[5].as<string>(),
                                 row[6].as<string>(),
                                 row[7].as<int64_t>()});

        return observations;
    };

    Observations SBSearchDatabasePostgreSQL::find_observations(const string &source, const double mjd_start, double mjd_stop, const int64_t limit, const int64_t offset)
    {
        error_if_closed();
        pqxx::nontransaction work(connection_);

        pqxx::result result = work.exec_params(
            "SELECT observation_id, source, observatory, product_id, mjd_start, mjd_stop, fov, terms "
            "FROM observations "
            "WHERE source = $1 AND mjd_start >= $2 AND mjd_stop <= $3 "
            "LIMIT $4 OFFSET $5",
            source,
            mjd_start,
            mjd_stop,
            limit,
            offset);

        Observations observations;
        observations.data.reserve(limit);
        for (auto const &row : result)
            observations.append({row[0].as<string>(),
                                 row[1].as<string>(),
                                 row[2].as<string>(),
                                 row[3].as<double>(),
                                 row[4].as<double>(),
                                 row[5].as<string>(),
                                 row[6].as<string>(),
                                 row[7].as<int64_t>()});

        return observations;
    };

    Observations SBSearchDatabasePostgreSQL::find_observations(vector<string> query_terms, const Options &options)
    {
        // query_terms may be spatial-temporal, just spatial, or just temporal.
        error_if_closed();
        pqxx::nontransaction work(connection_);

        int count = 0;
        std::set<int64_t> approximate_matches;

        // Query database with terms, but not too many at once
        vector<string> subset;
        subset.reserve(MAXIMUM_QUERY_TERMS);

        string statement = "SELECT observation_id FROM observations WHERE terms && $1";

        int index = 1;
        pqxx::params parameters;
        if (!options.source.empty())
        {
            statement.append(" AND source = $" + std::to_string(++index));
            parameters.append(options.source);
        }

        statement.append(" AND mjd_start >= $" + std::to_string(++index));
        parameters.append(options.mjd_start);

        statement.append(" AND mjd_stop <= $" + std::to_string(++index));
        parameters.append(options.mjd_stop);

        connection_.prepare("", statement);
        for (int i = 0; i < query_terms.size(); i += MAXIMUM_QUERY_TERMS)
        {
            subset.clear();

            const int j = std::min(query_terms.size(), i + MAXIMUM_QUERY_TERMS);
            std::copy(query_terms.begin() + i, query_terms.begin() + j, std::back_inserter(subset));

            pqxx::result result = work.exec_prepared("", subset, parameters);
            for (auto const &row : result)
                approximate_matches.insert(row[0].as<int64_t>());

            Logger::debug() << "Searched " << j << " of "
                            << query_terms.size() << " query terms."
                            << endl;
        }

        work.commit();

        return get_observations(approximate_matches.begin(), approximate_matches.end());
    };

    void SBSearchDatabasePostgreSQL::add_found(const Founds &founds)
    {
        Logger::info() << "Adding " << founds.size() << " found observations." << endl;
        error_if_closed();

        char now[32];
        std::time_t time_now = std::time(nullptr);
        std::strftime(now, 32, "%F %T", std::gmtime(&time_now));

        pqxx::work work(connection_);

        connection_.prepare("", R"(
            INSERT INTO found (
                observation_id,
                moving_target_id,
                mjd,
                tmtp,
                ra,
                dec,
                unc_a,
                unc_b,
                unc_theta,
                rh,
                delta,
                phase,
                selong,
                true_anomaly,
                sangle,
                vangle,
                vmag,
                saved
            ) VALUES (
                $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18
            ))");

        for (auto const &found : founds)
            work.exec_prepared(
                "",
                found.observation.observation_id(),
                found.ephemeris.target().moving_target_id().value(),
                found.ephemeris.data(0).mjd,
                found.ephemeris.data(0).tmtp,
                found.ephemeris.data(0).ra,
                found.ephemeris.data(0).dec,
                found.ephemeris.data(0).unc_a,
                found.ephemeris.data(0).unc_b,
                found.ephemeris.data(0).unc_theta,
                found.ephemeris.data(0).rh,
                found.ephemeris.data(0).delta,
                found.ephemeris.data(0).phase,
                found.ephemeris.data(0).selong,
                found.ephemeris.data(0).true_anomaly,
                found.ephemeris.data(0).sangle,
                found.ephemeris.data(0).vangle,
                found.ephemeris.data(0).vmag,
                now);

        work.commit();
    }

    Founds SBSearchDatabasePostgreSQL::get_found(const Observation &observation)
    {
        error_if_closed();

        pqxx::nontransaction work(connection_);
        pqxx::result result = work.exec_params(
            R"(
            SELECT
                moving_target_id, mjd, tmtp,
                ra, dec, unc_a, unc_b, unc_theta,
                rh, delta, phase, selong, true_anomaly,
                sangle, vangle, vmag
            FROM found
            WHERE observation_id=$1)",
            observation.observation_id());
        work.commit();

        Founds founds;
        for (auto const &row : result)
        {
            MovingTarget target = get_moving_target(row[0].as<int64_t>());

            Ephemeris::Datum d;
            d.mjd = row[1].as<double>();
            d.tmtp = row[2].as<double>();
            d.ra = row[3].as<double>();
            d.dec = row[4].as<double>();
            d.unc_a = row[5].as<double>();
            d.unc_b = row[6].as<double>();
            d.unc_theta = row[7].as<double>();
            d.rh = row[8].as<double>();
            d.delta = row[9].as<double>();
            d.phase = row[10].as<double>();
            d.selong = row[11].as<double>();
            d.true_anomaly = row[12].as<double>();
            d.sangle = row[13].as<double>();
            d.vangle = row[14].as<double>();
            d.vmag = row[15].as<double>();

            Ephemeris ephemeris(target, {d});
            founds.append(Found(observation, ephemeris));
        }
        return founds;
    };

    Founds SBSearchDatabasePostgreSQL::get_found(const MovingTarget &target)
    {
        error_if_closed();
        pqxx::nontransaction work(connection_);
        pqxx::result result = work.exec_params(
            R"(
            SELECT
                observation_id, mjd, tmtp,
                ra, dec, unc_a, unc_b, unc_theta,
                rh, delta, phase, selong, true_anomaly,
                sangle, vangle, vmag
            FROM found
            WHERE moving_target_id=$1)",
            target.moving_target_id());
        work.commit();

        Founds founds;
        for (auto const &row : result)
        {
            Observation observation = get_observation(row[0].as<int64_t>());

            Ephemeris::Datum d;
            d.mjd = row[1].as<double>();
            d.tmtp = row[2].as<double>();
            d.ra = row[3].as<double>();
            d.dec = row[4].as<double>();
            d.unc_a = row[5].as<double>();
            d.unc_b = row[6].as<double>();
            d.unc_theta = row[7].as<double>();
            d.rh = row[8].as<double>();
            d.delta = row[9].as<double>();
            d.phase = row[10].as<double>();
            d.selong = row[11].as<double>();
            d.true_anomaly = row[12].as<double>();
            d.sangle = row[13].as<double>();
            d.vangle = row[14].as<double>();
            d.vmag = row[15].as<double>();

            Ephemeris ephemeris(target, {d});
            founds.append(Found(observation, ephemeris));
        }

        return founds;
    };

    void SBSearchDatabasePostgreSQL::remove_found(const Founds &founds)
    {
        Logger::info() << "Removing " << founds.size() << " found entries." << endl;
        error_if_closed();

        pqxx::work work(connection_);
        connection_.prepare("", "DELETE FROM found WHERE moving_target_id=$1 and observation_id=$2");

        for (auto const found : founds)
        {
            if (!found.ephemeris.target().moving_target_id())
                throw MovingTargetError("Cannot remove found items for a moving target without an ID.");

            work.exec_prepared(
                "",
                found.ephemeris.target().moving_target_id(),
                found.observation.observation_id());
        }
        work.commit();
    };

    void SBSearchDatabasePostgreSQL::error_if_closed()
    {
        if (!connection_.is_open())
            throw std::runtime_error("Database is not open.");
    }
}