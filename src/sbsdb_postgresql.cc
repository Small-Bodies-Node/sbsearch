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
    SBSearchDatabasePostgreSQL::SBSearchDatabasePostgreSQL(const string uri)
    {
        connection_ = new pqxx::connection(uri);

        if (connection_->is_open())
            Logger::info() << "Opened database URI " << uri << endl;
        else
            throw "Failed to open database connection.";
    }

    void SBSearchDatabasePostgreSQL::close()
    {
        if (connection_->is_open())
        {
            connection_->close();
            Logger::info() << "Closed database." << endl;
        }
    };

    void SBSearchDatabasePostgreSQL::execute_sql(const char *statement)
    {
        error_if_closed();

        try
        {
            pqxx::nontransaction work(*connection_);
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

        pqxx::nontransaction work(*connection_);

        work.exec(R"(
CREATE TABLE IF NOT EXISTS observations (
  observation_id INTEGER PRIMARY KEY,
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

        work.exec(R"(
CREATE TABLE IF NOT EXISTS moving_targets (
  moving_targets_row_id INTEGER PRIMARY KEY,
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
  ephemeris_id INTEGER PRIMARY KEY,
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
  found_id INTEGER PRIMARY KEY,
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
  vmag DOUBLE PRECISION NOT NULL,
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

        Logger::debug() << "Database tables are set." << endl;
    };

    void SBSearchDatabasePostgreSQL::drop_observations_indices() {};
    void SBSearchDatabasePostgreSQL::create_observations_indices()
    {
        Logger::info() << "Creating observations indices." << std::endl;

        execute_sql(R"(
        CREATE INDEX IF NOT EXISTS ix_observations_terms
        ON observations
        USING GIN (terms);

        ANALYZE observations;
)");

        Logger::info()
            << "Created observations term index." << std::endl;
    };

    // get single value results from a SQL statement
    double *SBSearchDatabasePostgreSQL::get_double(const char *statement)
    {
        error_if_closed();

        pqxx::nontransaction work(*connection_);
        pqxx::row row = work.exec1(statement);
        double *value = new double;
        if (row[0].is_null())
            value = nullptr;
        else
            *value = row[0].as<double>();
        return std::move(value);
    };

    int *SBSearchDatabasePostgreSQL::get_int(const char *statement)
    {
        error_if_closed();
        pqxx::nontransaction work(*connection_);
        pqxx::row row = work.exec1(statement);
        int *value = new int;
        if (row[0].is_null())
            value = nullptr;
        else
            *value = row[0].as<int>();
        return std::move(value);
    };

    int64 *SBSearchDatabasePostgreSQL::get_int64(const char *statement)
    {
        error_if_closed();
        pqxx::nontransaction work(*connection_);
        pqxx::row row = work.exec1(statement);
        int64 *value = new int64;
        if (row[0].is_null())
            value = nullptr;
        else
            *value = row[0].as<int64>();
        return std::move(value);
    };

    string *SBSearchDatabasePostgreSQL::get_string(const char *statement)
    {
        error_if_closed();
        pqxx::nontransaction work(*connection_);
        pqxx::row row = work.exec1(statement);
        string *value = new string();
        if (row[0].is_null())
            value = nullptr;
        else
            *value = row[0].as<string>();
        return std::move(value);
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

        pqxx::work work(*connection_);
        for (int i = 0; i < parameters.size(); i++)
        {
            connection_->prepare("",
                                 "UPDATE configuration SET value = $1 WHERE parameter = $2");
            work.exec_prepared("", values[i], parameters[i]);
        }
        work.commit();
    };

    std::pair<double *, double *> SBSearchDatabasePostgreSQL::observation_date_range(const string &source) {};

    void SBSearchDatabasePostgreSQL::add_moving_target(MovingTarget &target) {};
    void SBSearchDatabasePostgreSQL::remove_moving_target(const MovingTarget &target) {};
    void SBSearchDatabasePostgreSQL::update_moving_target(const MovingTarget &target) {};
    MovingTarget SBSearchDatabasePostgreSQL::get_moving_target(const int moving_target_id) {};
    MovingTarget SBSearchDatabasePostgreSQL::get_moving_target(const string &name, const bool small_body) {};
    vector<MovingTarget> SBSearchDatabasePostgreSQL::get_all_moving_targets() {};

    void SBSearchDatabasePostgreSQL::add_observatory(const string &name, const Observatory &observatory) {};
    const Observatory SBSearchDatabasePostgreSQL::get_observatory(const string &name) {};
    const Observatories SBSearchDatabasePostgreSQL::get_observatories() {};
    void SBSearchDatabasePostgreSQL::remove_observatory(const string &name) {};
    const vector<string> SBSearchDatabasePostgreSQL::get_sources() {};

    void SBSearchDatabasePostgreSQL::add_ephemeris(Ephemeris &eph) {};
    Ephemeris SBSearchDatabasePostgreSQL::get_ephemeris(const MovingTarget target, double mjd_start, double mjd_stop) {};
    int SBSearchDatabasePostgreSQL::remove_ephemeris(const MovingTarget target, double mjd_start, double mjd_stop) {};
    std::pair<double *, double *> SBSearchDatabasePostgreSQL::ephemeris_date_range() {};

    void SBSearchDatabasePostgreSQL::add_observation(Observation &observation) {};
    Observation SBSearchDatabasePostgreSQL::get_observation(const int64 observation_id) {};
    void SBSearchDatabasePostgreSQL::remove_observations(const double mjd_start, const double mjd_stop) {};
    void SBSearchDatabasePostgreSQL::remove_observations(const string &source, const double mjd_start, const double mjd_stop) {};

    // Count number of observations within an interval.
    int64 SBSearchDatabasePostgreSQL::count_observations(const double mjd_start, const double mjd_stop) {};

    // Count number of observations for a source within an interval, if
    // source is an empty string, then count all sources.
    int64 SBSearchDatabasePostgreSQL::count_observations(const string &source, const double mjd_start, const double mjd_stop) {};

    Observations SBSearchDatabasePostgreSQL::find_observations(const double mjd_start, const double mjd_stop, const int64 limit, const int64 offset) {};
    Observations SBSearchDatabasePostgreSQL::find_observations(const string &source, const double mjd_start, double mjd_stop, const int64 limit, const int64 offset) {};
    Observations SBSearchDatabasePostgreSQL::find_observations(vector<string> query_terms, const Options &options) {};

    void SBSearchDatabasePostgreSQL::add_found(const Found &found) {};
    Founds SBSearchDatabasePostgreSQL::get_found(const Observation &observation) {};
    Founds SBSearchDatabasePostgreSQL::get_found(const MovingTarget &target) {};
    void SBSearchDatabasePostgreSQL::remove_found(const Found &found) {};

    void SBSearchDatabasePostgreSQL::error_if_closed()
    {
        if (!connection_->is_open())
            throw std::runtime_error("Database is not open.");
    }
}