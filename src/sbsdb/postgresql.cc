#include "config.h"

#include <optional>
#include <pqxx/pqxx>

#include "./postgresql.h"
#include "ephemeris.h"
#include "exceptions.h"
#include "logging.h"
#include "observation.h"

using std::endl;

namespace sbsearch::sbsdb
{
    bool Postgresql::begin()
    {
        if (in_transaction_)
            return false;

        execute("BEGIN");
        in_transaction_ = true;
        return true;
    }

    void Postgresql::rollback()
    {
        if (!in_transaction_)
            return;

        execute("ROLLBACK");
        in_transaction_ = false;
    }

    void Postgresql::commit()
    {
        if (!in_transaction_)
            return;

        execute("COMMIT");
        in_transaction_ = false;
    }

    template <>
    int Postgresql::row_as(const pqxx::row &row, const int i)
    {
        return row[i].as<int>();
    }

    template <>
    int64_t Postgresql::row_as(const pqxx::row &row, const int i)
    {
        return row[i].as<int64_t>();
    }

    template <>
    double Postgresql::row_as(const pqxx::row &row, const int i)
    {
        return row[i].as<double>();
    }

    template <>
    string Postgresql::row_as(const pqxx::row &row, const int i)
    {
        return row[i].as<string>();
    }

    template <>
    optional<int> Postgresql::row_as(const pqxx::row &row, const int i)
    {
        if (row[i].is_null())
            return {};

        return {row[i].as<int>()};
    }

    template <>
    optional<int64_t> Postgresql::row_as(const pqxx::row &row, const int i)
    {
        if (row[i].is_null())
            return {};

        return {row[i].as<int64_t>()};
    }

    template <>
    optional<double> Postgresql::row_as(const pqxx::row &row, const int i)
    {
        if (row[i].is_null())
            return {};

        return {row[i].as<double>()};
    }

    template <>
    optional<string> Postgresql::row_as(const pqxx::row &row, const int i)
    {
        if (row[i].is_null())
            return {};

        return {row[i].as<string>()};
    }

    template <>
    Ephemeris::Datum Postgresql::row_as(const pqxx::row &row, const int i)
    {
        Ephemeris::Datum d;
        d.mjd = row[row.column_number("mjd")].as<double>();
        d.tmtp = row[row.column_number("tmtp")].as<double>();
        d.ra = row[row.column_number("ra")].as<double>();
        d.dec = row[row.column_number("dec")].as<double>();
        d.mu = row[row.column_number("mu")].as<double>();
        d.mu_theta = row[row.column_number("mu_theta")].as<double>();
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
    Found::DBModel Postgresql::row_as(const pqxx::row &row, const int i)
    {
        Found::DBModel model;
        model.found_id = row[row.column_number("found_id")].as<int64_t>();
        model.observation_id = row[row.column_number("observation_id")].as<int64_t>();
        model.moving_target_id = row[row.column_number("moving_target_id")].as<int64_t>();
        model.mjd = row[row.column_number("mjd")].as<double>();
        model.tmtp = row[row.column_number("tmtp")].as<double>();
        model.ra = row[row.column_number("ra")].as<double>();
        model.dec = row[row.column_number("dec")].as<double>();
        model.mu = row[row.column_number("mu")].as<double>();
        model.mu_theta = row[row.column_number("mu_theta")].as<double>();
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
        model.mjd_added = row[row.column_number("mjd_added")].as<double>();
        return model;
    }

    template <>
    MovingTarget::DBModel Postgresql::row_as(const pqxx::row &row, const int i)
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
    Observation Postgresql::row_as(const pqxx::row &row, const int i)
    {
        const int64_t observation_id = row[row.column_number("observation_id")].as<int64_t>();
        const string source = row[row.column_number("source")].as<string>();
        const string observatory = row[row.column_number("observatory")].as<string>();
        const string product_id = row[row.column_number("product_id")].as<string>();
        const double mjd_start = row[row.column_number("mjd_start")].as<double>();
        const double mjd_stop = row[row.column_number("mjd_stop")].as<double>();
        const string fov = row[row.column_number("fov")].as<string>();
        const string center = row[row.column_number("center")].as<string>();
        const optional<string> meta = row[row.column_number("meta")].as<optional<string>>();
        const double mjd_added = row[row.column_number("mjd_added")].as<double>();

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

        return {
            source,
            observatory,
            product_id,
            mjd_start,
            mjd_stop,
            fov,
            terms,
            observation_id,
            center,
            meta,
            mjd_added,
        };
    }

    template <>
    Observatory Postgresql::row_as(const pqxx::row &row, const int i)
    {
        const string name = row[row.column_number("name")].as<string>();
        const double longitude = row[row.column_number("longitude")].as<double>();
        const double rho_cos_phi = row[row.column_number("rho_cos_phi")].as<double>();
        const double rho_sin_phi = row[row.column_number("rho_sin_phi")].as<double>();

        return {longitude, rho_cos_phi, rho_sin_phi, name};
    }

    vector<Observation> Postgresql::get_all_observations(const string &table)
    {
        const int count = get_one<int>("SELECT COUNT(*) FROM " + work_.quote_name(table));
        if (count == 0)
            return {};

        vector<Observation> observations;
        observations.reserve(count);

        auto stream = work_.stream<string, string, string, double, double, string, vector<string>, int64_t, string, optional<string>, double>(
            "SELECT DISTINCT ON (observation_id) source,observatory,product_id,mjd_start,mjd_stop,"
            "fov,terms,observation_id,center,meta,mjd_added FROM " +
            work_.quote_name(table));

        for (auto row : stream)
            std::apply([&](auto... args)
                       { observations.emplace_back(args...); }, row);

        return observations;
    }

    void Postgresql::insert_ephemeris(const Ephemeris &ephemeris)
    {
        auto insert = pqxx::stream_to::table(
            work_,
            {"ephemerides"},
            {"moving_target_id", "mjd", "tmtp", "ra", "dec",
             "mu", "mu_theta", "unc_a", "unc_b", "unc_theta",
             "rh", "delta", "phase", "selong", "true_anomaly",
             "sangle", "vangle", "vmag", "mjd_added"});

        for (auto const &row : ephemeris.data())
        {
            insert.write_values(ephemeris.target().moving_target_id(),
                                row.mjd, row.tmtp, row.ra, row.dec,
                                row.mu, row.mu_theta, row.unc_a, row.unc_b, row.unc_theta,
                                row.rh, row.delta, row.phase, row.selong, row.true_anomaly,
                                row.sangle, row.vangle, row.vmag, Date::now().mjd());
        }
        insert.complete();
    }

    size_t Postgresql::insert_many_observations(Observations &observations)
    {
        // insert into a temporary table, then we insert those results into
        // the observation table, returning the new observation_ids.
        execute("CREATE TEMPORARY TABLE insert_observations (LIKE observations)");
        execute("ALTER TABLE insert_observations DROP COLUMN observation_id");

        auto insert = pqxx::stream_to::table(
            work_,
            {"insert_observations"},
            {"source", "observatory", "product_id", "mjd_start", "mjd_stop",
             "fov", "terms", "center", "meta", "mjd_added"});

        for (auto const &obs : observations)
        {
            string terms = "{" + util::join(obs.terms(), ",") + "}";
            insert.write_values(obs.source(), obs.observatory(), obs.product_id(),
                                obs.mjd_start(), obs.mjd_stop(), obs.fov(), terms,
                                obs.center(), obs.meta(), obs.mjd_added());
        }
        insert.complete();

        // update observations table, returning observation_id
        auto returning = work_.stream<string, int64_t>(
            "INSERT INTO observations "
            "(source, observatory, product_id, mjd_start, mjd_stop, fov, terms, center, meta, mjd_added) "
            "(SELECT source, observatory, product_id, mjd_start, mjd_stop, fov, terms, center, meta, mjd_added "
            " FROM insert_observations) "
            "RETURNING product_id, observation_id");

        std::map<string, int64_t> observation_id;
        for (auto const &[pid, oid] : returning)
            observation_id[pid] = oid;

        for (auto &obs : observations)
            obs.observation_id(observation_id[obs.product_id()]);

        execute("DROP TABLE insert_observations");

        return observation_id.size();
    }

    size_t Postgresql::update_many_observation_terms(const vector<int64_t> &observation_ids,
                                                     const vector<vector<string>> &observation_terms)
    {
        assert(observation_ids.size() == observation_terms.size());

        // insert into a temporary table, then we insert those results into
        // the observation table.
        execute("CREATE TEMPORARY TABLE update_observation_terms "
                "(observation_id INTEGER PRIMARY KEY,"
                "terms TEXT[] NOT NULL)");

        auto update = pqxx::stream_to::table(
            work_,
            {"update_observation_terms"},
            {"observation_id", "terms"});

        for (size_t i = 0; i < observation_ids.size(); i++)
        {
            string terms = "{" + util::join(observation_terms[i], ",") + "}";
            update.write_values(observation_ids[i], terms);
        }
        update.complete();

        // update observations table, returning observation_id
        execute("UPDATE observations AS a SET terms = b.terms "
                "FROM update_observation_terms AS b "
                "WHERE a.observation_id = b.observation_id");

        execute("DROP TABLE update_observation_terms");

        return observation_ids.size();
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
  fov VARCHAR(1024) NOT NULL,
  center VARCHAR(16) NOT NULL,
  terms TEXT[] NOT NULL,
  meta TEXT[],
  mjd_added DOUBLE PRECISION NOT NULL
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
  mu DOUBLE PRECISION NOT NULL,
  mu_theta DOUBLE PRECISION NOT NULL,
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
  mjd_added DOUBLE PRECISION NOT NULL
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
  mu DOUBLE PRECISION NOT NULL,
  mu_theta DOUBLE PRECISION NOT NULL,
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
  mjd_added DOUBLE PRECISION NOT NULL
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
INSERT INTO configuration VALUES ('database version', ')" SBSEARCH_VERSION R"(') ON CONFLICT DO NOTHING;

ANALYZE;
)");

        Logger::debug() << "Database tables are set." << endl;
    }

    void Postgresql::drop_observations_indices()
    {
        Logger::info() << "Dropping observations indices." << std::endl;
        execute("DROP INDEX IF EXISTS idx_observations_terms;");
        execute("DROP INDEX IF EXISTS idx_observations_mjd_start;");
        execute("DROP INDEX IF EXISTS idx_observations_mjd_stop;");
        execute("DROP INDEX IF EXISTS idx_observations_source_mjd_start;");
        execute("DROP INDEX IF EXISTS idx_observations_source_mjd_stop;");
        execute("DROP INDEX IF EXISTS idx_observations_observatory;");
        execute("DROP INDEX IF EXISTS idx_observations_product_id;");
        execute("DROP INDEX IF EXISTS idx_observations_center;");
        Logger::info() << "Observations indices dropped." << std::endl;
    };

    void Postgresql::create_observations_indices()
    {
        Logger::info() << "Creating observations table indices." << std::endl;

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

        CREATE INDEX IF NOT EXISTS idx_observations_center
        ON observations(center);

        CLUSTER observations USING idx_observations_center;

        ANALYZE observations;
)");

        Logger::info()
            << "Created observations indices." << std::endl;
    };
};