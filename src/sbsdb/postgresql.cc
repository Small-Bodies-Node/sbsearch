#include <optional>
#include <pqxx/pqxx>

#include "sbsdb.h"
#include "../ephemeris.h"
#include "../exceptions.h"
#include "../observation.h"
#include "postgresql.h"

namespace sbsearch::sbsdb
{
    template <>
    int Postgresql::row_as(const pqxx::row &row)
    {
        return row[0].as<int>();
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
    optional<string> Postgresql::row_as(const pqxx::row &row)
    {
        if (row[0].is_null())
            return {};

        return {row[0].as<string>()};
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
};