#ifndef SBSDB_POSTGRESQL_H_
#define SBSDB_POSTGRESQL_H_

#include <optional>
#include <string>
#include <pqxx/pqxx>

#include "../observation.h"
#include "../util.h"

using std::optional;
using std::string;

namespace sbsearch::sbsdb
{
    class Postgresql
    {
    public:
        Postgresql(const string &url) : connection_(url) {};

        // Execute a statement, parameters are optional, nothing is returned.
        template <typename... Targs>
        void execute(const string &statement, Targs... args)
        {
            const int nargs = sizeof...(args);

            pqxx::transaction work(connection_);

            if (nargs == 0)
                work.exec(statement);
            else
                work.exec_params(statement, pqxx::params(args...));
        }

        // Execute a statement, parameters are optional, return a single value
        // from the first column of the first row.
        template <typename T, typename... Targs>
        T get_one(const string &statement, Targs... args)
        {
            const int nargs = sizeof...(args);

            pqxx::transaction work(connection_);
            pqxx::row row;

            if (nargs == 0)
                row = work.exec1(statement);
            else
                row = work.exec_params1(statement, pqxx::params(args...));
            return row_as<T>(row);
        }

        // Execute a statement, parameters are optional, return a vector of values
        template <typename T, typename... Targs>
        auto get_many(const string &statement, Targs... args)
            -> std::enable_if_t<std::is_arithmetic_v<T>, std::vector<T>>
        {
            const int nargs = sizeof...(args);

            pqxx::transaction work(connection_);
            pqxx::result result;

            if (nargs == 0)
                result = work.exec(statement);
            else
                result = work.exec_params(statement, pqxx::params(args...));

            std::vector<T> v(result.size());
            std::transform(result.begin(), result.end(), v.begin(),
                           [&](const pqxx::row &row)
                           { return row[0].as<T>(); });

            return v;
        }

        // Execute a statement, parameters are optional, returning an
        // Observations object.
        template <typename T, typename... Targs>
        auto get_many(const string &statement, Targs... args)
            -> std::enable_if_t<std::is_same_v<T, Observation>, Observations>
        {
            const int nargs = sizeof...(args);

            pqxx::transaction work(connection_);
            pqxx::result result;

            if (nargs == 0)
                result = work.exec(statement);
            else
                result = work.exec_params(statement, pqxx::params(args...));

            Observations observations;
            observations.data.resize(result.size());
            std::transform(result.begin(), result.end(), observations.data.begin(),
                           [&](const pqxx::row &row)
                           { return row_as<Observation>(row); });

            return observations;
        }

        // Execute a statement, parameters are optional, returning an
        // Ephemeris object.
        template <typename T, typename... Targs>
        auto get_many(const string &statement, Targs... args)
            -> std::enable_if_t<std::is_same_v<T, Ephemeris::Datum>, Ephemeris::Data>
        {
            const int nargs = sizeof...(args);

            pqxx::transaction work(connection_);
            pqxx::result result;

            if (nargs == 0)
                result = work.exec(statement);
            else
                result = work.exec_params(statement, pqxx::params(args...));

            Ephemeris::Data data(result.size());
            std::transform(result.begin(), result.end(), data.begin(),
                           [&](const pqxx::row &row)
                           { return row_as<Ephemeris::Datum>(row); });

            return data;
        }

    private:
        pqxx::connection connection_;

        template <typename T>
        T row_as(const pqxx::row &row);
    };
}

#endif // SBSDB_POSTGRESQL_H_
