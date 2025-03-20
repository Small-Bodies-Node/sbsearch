#ifndef SBSDB_POSTGRESQL_H_
#define SBSDB_POSTGRESQL_H_

#include <optional>
#include <string>
#include <pqxx/pqxx>

#include "sbsdb.h"
#include "../ephemeris.h"
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

        /**
         * @brief Execute an SQL statement.
         *
         * @param statement The statement to execute.
         * @param args Optional parameters.
         */
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

        /**
         * @brief Execute an SQL statement returning a single value.
         *
         * @tparam T The data type.  May be std::optional<int>,
         *           std::optional<int64_t>, std::optional<string>, Observation,
         *           Ephemeris::Datum, or MovingTargetsModel.
         * @param statement The statement to execute.  The "terms" column of
         *                  "observations" is optional for the Observation type.
         * @param args Optional parameters.
         * @return T The retrieved value.  Basic types are derived from the
         *           first column, sbsearch types from the full row.
         */
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

        /**
         * @brief Execute an SQL statement returning a vector of values.
         *
         * @tparam T The vector data type.
         * @param statement The statement to execute.
         * @param args Optional parameters.
         * @return std::vector<T> The values from the first column of the results.
         */
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

        /**
         * @brief Execute an SQL statement returning an Observations object.
         *
         * @tparam T Observation
         * @param statement The statement to execute.
         * @param args Optional parameters.
         * @return Observations The observations based on each row of the results.
         */
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

        /**
         * @brief Execute an SQL statement returning an Ephemeris::Data object.
         *
         * @tparam T Ephemeris::Data
         * @param statement The statement to execute.
         * @param args Optional parameters.
         * @return Ephemeris::Data The ephemeris data based on each row of the results.
         */
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

        /**
         * @brief Execute an SQL statement returning a vector of moving target
         * data.
         *
         * @tparam T MovingTargetsModel
         * @param statement The statement to execute.
         * @param args Optional parameters.
         * @return std::vector<MovingTargetsModel> The moving target data based
         *                                         on each row of the results.
         */
        template <typename T, typename... Targs>
        auto get_many(const string &statement, Targs... args)
            -> std::enable_if_t<std::is_same_v<T, MovingTargetsModel>, std::vector<MovingTargetsModel>>
        {
            const int nargs = sizeof...(args);

            pqxx::transaction work(connection_);
            pqxx::result result;

            if (nargs == 0)
                result = work.exec(statement);
            else
                result = work.exec_params(statement, pqxx::params(args...));

            std::vector<MovingTargetsModel> rows(result.size());
            std::transform(result.begin(), result.end(), rows.begin(),
                           [&](const pqxx::row &row)
                           { return row_as<MovingTargetsModel>(row); });

            return rows;
        }

    private:
        pqxx::connection connection_;

        template <typename T>
        T row_as(const pqxx::row &row);
    };
}

#endif // SBSDB_POSTGRESQL_H_
