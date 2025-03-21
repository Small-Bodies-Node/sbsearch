#include <optional>
#include <pqxx/pqxx>

#include "sbsdb.h"
#include "postgresql.h"
#include "../config.h"
#include "../ephemeris.h"
#include "../exceptions.h"
#include "../logging.h"
#include "../observation.h"

using std::endl;

namespace sbsearch::sbsdb
{
    void Postgresql::begin()
    {
        // Logger::debug() << "Begin database transaction." << endl;
        execute("BEGIN");
        in_transaction_ = true;
    }

    void Postgresql::rollback()
    {
        // Logger::debug() << "Rollback database transaction." << endl;
        execute("ROLLBACK");
        in_transaction_ = false;
    }

    void Postgresql::commit()
    {
        // Logger::debug() << "Commit database transaction." << endl;
        execute("COMMIT");
        in_transaction_ = false;
    }

    template <>
    int Postgresql::row_as(const pqxx::row &row)
    {
        return row[0].as<int>();
    }

    template <>
    int64_t Postgresql::row_as(const pqxx::row &row)
    {
        return row[0].as<int64_t>();
    }

    template <>
    double Postgresql::row_as(const pqxx::row &row)
    {
        return row[0].as<double>();
    }

    template <>
    string Postgresql::row_as(const pqxx::row &row)
    {
        return row[0].as<string>();
    }

    template <>
    optional<int> Postgresql::row_as(const pqxx::row &row)
    {
        if (row[0].is_null())
            return {};

        return {row[0].as<int>()};
    }

    template <>
    optional<int64_t> Postgresql::row_as(const pqxx::row &row)
    {
        if (row[0].is_null())
            return {};

        return {row[0].as<int64_t>()};
    }

    template <>
    optional<double> Postgresql::row_as(const pqxx::row &row)
    {
        if (row[0].is_null())
            return {};

        return {row[0].as<double>()};
    }

    template <>
    optional<string> Postgresql::row_as(const pqxx::row &row)
    {
        if (row[0].is_null())
            return {};

        return {row[0].as<string>()};
    }

    template <>
    Ephemeris::Datum Postgresql::row_as(const pqxx::row &row)
    {
        Ephemeris::Datum d;
        d.mjd = row[row.column_number("mjd")].as<double>();
        d.tmtp = row[row.column_number("tmtp")].as<double>();
        d.ra = row[row.column_number("ra")].as<double>();
        d.dec = row[row.column_number("dec")].as<double>();
        d.unc_a = row[row.column_number("unc_a")].as<double>();
        d.unc_b = row[row.column_number("unc_b")].as<double>();
        d.unc_theta = row[row.column_number("unc_theta")].as<double>();
        d.rh = row[row.column_number("rh")].as<double>();
        d.delta = row[row.column_number("delta")].as<double>();
        d.phase = row[row.column_number("phase")].as<double>();
        d.selong = row[row.column_number("selong")].as<double>();
        d.true_anomaly = row[row.column_number("true_anomaly")].as<double>();
        d.sangle = row[row.column_number("sangle")].as<double>();
        d.vangle = row[row.column_number("vangle")].as<double>();
        d.vmag = row[row.column_number("vmag")].as<double>();
        return d;
    }

    template <>
    Found::DBModel Postgresql::row_as(const pqxx::row &row)
    {
        Found::DBModel model;
        model.found_id = row[row.column_number("found_id")].as<int64_t>();
        model.observation_id = row[row.column_number("observation_id")].as<int64_t>();
        model.moving_target_id = row[row.column_number("moving_target_id")].as<int64_t>();
        model.mjd = row[row.column_number("mjd")].as<double>();
        model.tmtp = row[row.column_number("tmtp")].as<double>();
        model.ra = row[row.column_number("ra")].as<double>();
        model.dec = row[row.column_number("dec")].as<double>();
        model.unc_a = row[row.column_number("unc_a")].as<double>();
        model.unc_b = row[row.column_number("unc_b")].as<double>();
        model.unc_theta = row[row.column_number("unc_theta")].as<double>();
        model.rh = row[row.column_number("rh")].as<double>();
        model.delta = row[row.column_number("delta")].as<double>();
        model.phase = row[row.column_number("phase")].as<double>();
        model.selong = row[row.column_number("selong")].as<double>();
        model.true_anomaly = row[row.column_number("true_anomaly")].as<double>();
        model.sangle = row[row.column_number("sangle")].as<double>();
        model.vangle = row[row.column_number("vangle")].as<double>();
        model.vmag = row[row.column_number("vmag")].as<double>();
        model.saved = row[row.column_number("saved")].as<string>();
        return model;
    }

    template <>
    MovingTarget::DBModel Postgresql::row_as(const pqxx::row &row)
    {
        MovingTarget::DBModel model;
        model.moving_targets_row_id = row[row.column_number("moving_targets_row_id")].as<int64_t>();
        model.moving_target_id = row[row.column_number("moving_target_id")].as<int64_t>();
        model.name = row[row.column_number("name")].as<string>();
        model.small_body = row[row.column_number("small_body")].as<bool>();
        model.primary_id = row[row.column_number("primary_id")].as<bool>();
        return model;
    }

    template <>
    Observation Postgresql::row_as(const pqxx::row &row)
    {
        const int64_t observation_id = row[row.column_number("observation_id")].as<int64_t>();
        const string source = row[row.column_number("source")].as<string>();
        const string observatory = row[row.column_number("observatory")].as<string>();
        const string product_id = row[row.column_number("product_id")].as<string>();
        const double mjd_start = row[row.column_number("mjd_start")].as<double>();
        const double mjd_stop = row[row.column_number("mjd_stop")].as<double>();
        const string fov = row[row.column_number("fov")].as<string>();

        // terms is optional
        vector<string> terms;
        try
        {
            int i = row.column_number("terms");
            auto parsed = row[i].as_array();

            std::pair<pqxx::array_parser::juncture, string> next;
            do
            {
                next = parsed.get_next();
                if (next.first == pqxx::array_parser::juncture::string_value)
                    terms.push_back(next.second);
            } while (next.first != pqxx::array_parser::juncture::done);
        }
        catch (pqxx::argument_error &e)
        {
        }

        return {source, observatory, product_id, mjd_start, mjd_stop, fov, terms, observation_id};
    }

    template <>
    Observatory Postgresql::row_as(const pqxx::row &row)
    {
        const string name = row[row.column_number("name")].as<string>();
        const double longitude = row[row.column_number("longitude")].as<double>();
        const double rho_cos_phi = row[row.column_number("rho_cos_phi")].as<double>();
        const double rho_sin_phi = row[row.column_number("rho_sin_phi")].as<double>();

        return {longitude, rho_cos_phi, rho_sin_phi, name};
    }

    void Postgresql::setup_tables()
    {
        execute(R"(
CREATE TABLE IF NOT EXISTS observations (
  observation_id INTEGER PRIMARY KEY GENERATED BY DEFAULT AS IDENTITY,
  source TEXT NOT NULL,
  observatory VARCHAR(64) NOT NULL,
  product_id VARCHAR(128) NOT NULL,
  mjd_start DOUBLE PRECISION NOT NULL,
  mjd_stop DOUBLE PRECISION NOT NULL,
  fov VARCHAR(128) NOT NULL,
  terms TEXT[] NOT NULL
);
)");

        create_observations_indices();

        execute(R"(
CREATE TABLE IF NOT EXISTS moving_targets (
  moving_targets_row_id INTEGER PRIMARY KEY GENERATED BY DEFAULT AS IDENTITY,
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
  ephemeris_id INTEGER PRIMARY KEY GENERATED BY DEFAULT AS IDENTITY,
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
  found_id INTEGER PRIMARY KEY GENERATED BY DEFAULT AS IDENTITY,
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
INSERT INTO configuration VALUES ('database version', ')" SBSEARCH_VERSION R"(') ON CONFLICT DO NOTHING;

ANALYZE;
)");

        Logger::debug() << "Database tables are set." << endl;
    }

    void Postgresql::drop_observations_indices()
    {
        Logger::info() << "Dropping observations indices." << std::endl;
        execute("DROP INDEX idx_observations_terms;");
        execute("DROP INDEX idx_observations_mjd_start;");
        execute("DROP INDEX idx_observations_mjd_stop;");
        execute("DROP INDEX idx_observations_source_mjd_start;");
        execute("DROP INDEX idx_observations_source_mjd_stop;");
        execute("DROP INDEX idx_observations_observatory;");
        execute("DROP INDEX idx_observations_product_id;");
        Logger::info() << "Observations indices dropped." << std::endl;
    };

    void Postgresql::create_observations_indices()
    {
        Logger::info() << "Creating observations indices." << std::endl;

        execute(R"(
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
};