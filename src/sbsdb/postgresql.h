#ifndef SBSDB_POSTGRESQL_H_
#define SBSDB_POSTGRESQL_H_

#include <iostream>
#include <optional>
#include <string>
#include <sstream>
#include <string_view>
#include <tuple>
#include <vector>
#include <pqxx/pqxx>

#include "sbsdb.h"
#include "../ephemeris.h"
#include "../exceptions.h"
#include "../found.h"
#include "../moving_target.h"
#include "../observation.h"
#include "../util/string.h"

using std::optional;
using std::string;
using std::vector;

namespace sbsearch::sbsdb
{
    class Postgresql
    {
    public:
        /**
         * @brief Construct a new Postgresql object
         *
         * @param uri The parameters of the database connection, e.g.,
         *            "postgres:///sbsearch".
         */
        Postgresql(const string &uri) : connection_(uri) {};

        /**
         * @brief Execute an SQL statement.
         *
         * @param statement The statement to execute.
         *
         * @param args Optional parameters.
         */
        template <typename... Targs>
        void execute(const string &statement, Targs... args)
        {
            const int nargs = sizeof...(args);

            if (nargs == 0)
                work_.exec(statement);
            else
                work_.exec_params(statement, pqxx::params(args...));
        }

        /**
         * @brief Begin a transaction.
         *
         * @return true If a new transaction is started.
         *
         * @return false If a transaction was already open.
         *
         */
        bool begin();

        /**
         * @brief Is the db in a transaction?
         *
         * @return true
         * @return false
         */
        inline bool in_transaction() { return in_transaction_; }

        /**
         * @brief Rollback the transaction.
         */
        void rollback();

        /**
         * @brief Commit the transaction.
         */
        void commit();

        /**
         * @brief Efficiently add data to the database.
         *
         * @param statement The SQL statement to execute.
         *
         * @param data The data to add.
         */
        // template <typename... Targs>
        // void add_many(const string &statement, Targs... args);

        /**
         * @brief Execute an SQL statement returning a single value.
         *
         * @tparam T The data type.  May be std::optional<int>,
         *           std::optional<int64_t>, std::optional<string>, Observation,
         *           Observatory, Ephemeris::Datum, or MovingTargetsModel.
         *
         * @param statement The statement to execute.  The "terms" column of
         *                  "observations" is optional for the Observation type.
         *                  Must only return one row.
         *
         * @param args Optional parameters.
         *
         * @return T The retrieved value.  Basic types are derived from the
         *           first column, sbsearch types from the full row.
         */
        template <typename T, typename... Targs>
        T get_one(const string &statement, Targs... args)
        {
            const int nargs = sizeof...(args);

            pqxx::row row;

            if (nargs == 0)
                row = work_.exec1(statement);
            else
                row = work_.exec_params1(statement, pqxx::params(args...));

            try
            {
                return row_as<T>(row);
            }
            catch (pqxx::argument_error &e)
            {
                throw SBSException(e.what());
            }
        }

        /**
         * @brief Execute an SQL statement returning a vector of values.
         *
         * The values are based on the first column of the first row for simple
         * data types like int, string, etc., or based on the entire row for
         * sbsearch data types like Observatory.
         *
         * @tparam T The vector data type.
         *
         * @param statement The statement to execute.
         *
         * @param args Optional parameters.
         *
         * @return vector<T>
         */
        template <typename T, typename... Targs>
        auto get_many(const string &statement, Targs... args)
            -> std::enable_if_t<std::is_arithmetic_v<T> |
                                    std::is_same_v<T, std::string> |
                                    std::is_same_v<T, Ephemeris::Datum> |
                                    std::is_same_v<T, Found::DBModel> |
                                    std::is_same_v<T, MovingTarget::DBModel> |
                                    std::is_same_v<T, Observation> |
                                    std::is_same_v<T, Observatory>,
                                vector<T>>
        {
            // const int nargs = sizeof...(args);

            pqxx::result result;

            // if (nargs == 0)
            //     result = work_.exec(statement);
            // else
            //     result = work_.exec_params(statement, pqxx::params(args...));
            result = work_.exec_params(statement, pqxx::params(args...));

            vector<T> v(result.size());
            std::transform(result.begin(), result.end(), v.begin(),
                           [&](const pqxx::row &row)
                           { return row_as<T>(row); });

            return v;
        }

        /**
         * @brief Get all observations from a table.
         */
        vector<Observation> get_all_observations(const string &table)
        {
            const int count = get_one<int>("SELECT COUNT(*) FROM " + work_.quote_name(table));
            if (count == 0)
                return {};

            vector<Observation> observations;
            observations.reserve(count);

            auto stream = work_.stream<string, string, string, double, double, string, vector<string>, int64_t, string>(
                "SELECT DISTINCT ON (observation_id) source,observatory,product_id,mjd_start,mjd_stop,"
                "fov,terms,observation_id,center FROM " +
                work_.quote_name(table));

            for (auto row : stream)
                std::apply([&](auto... args)
                           { observations.emplace_back(args...); }, row);

            return observations;
        }

        /**
         * @brief Execute an SQL statement returning a vector of values.
         *
         * The values returned are a pair are based on the first and second
         * columns, respectively.
         *
         * @tparam T The first column data type.
         *
         * @tparam U The second column data type.
         *
         * @param statement The statement to execute.
         *
         * @param args Optional parameters.
         *
         * @return vector<std::pair<T, U>> The retrieved values.
         */
        template <typename T, typename U, typename... Targs>
        vector<std::pair<T, U>> get_many(const string &statement, Targs... args)
        {
            const int nargs = sizeof...(args);

            pqxx::result result;

            if (nargs == 0)
                result = work_.exec(statement);
            else
                result = work_.exec_params(statement, pqxx::params(args...));

            vector<std::pair<T, U>> v(result.size());
            std::transform(result.begin(), result.end(), v.begin(),
                           [&](const pqxx::row &row)
                           { return std::make_pair(row_as<T>(row, 0), row_as<U>(row, 1)); });
            return v;
        }

        /**
         * @brief Set up database tables.
         *
         */
        void setup_tables();

        /**
         * @brief Drop observations indices.
         *
         */
        void drop_observations_indices();

        /**
         * @brief Create observations indices.
         *
         */
        void create_observations_indices();

    private:
        pqxx::connection connection_;
        pqxx::nontransaction work_{connection_};
        bool in_transaction_ = false;

        template <typename T>
        T row_as(const pqxx::row &row, const int i = 0);
    };
}

template <>
struct pqxx::string_traits<vector<string>>
{
    static constexpr bool converts_to_string{false};
    static constexpr bool converts_from_string{true};

    // Write string version into buffer, or return a constant string.
    // static zview to_buf(char *begin, char *end, T const &value);

    // Write string version into buffer.
    static char *into_buf(char *begin, char *end, vector<string> const &value)
    {
        return string_traits<string>::into_buf(begin, end, to_string(value));
    }

    // Converting value to string may require this much buffer space at most.
    static std::size_t size_buffer(vector<string> const &value) noexcept
    {
        auto add = [](int x, std::string_view s)
        { return x + s.length() + 1; };

        return std::accumulate(value.begin(), value.end(), 0, add) + 2;
    }

    static vector<string> from_string(std::string_view text)
    {
        // psql array format is "{a,b,c}"
        if ((text[0] != '{') | (text[text.length() - 1] != '}'))
            throw pqxx::conversion_error("text does not appear to be a PSQL array");

        return sbsearch::util::split(text.substr(1, text.length() - 2), ',');
    }

    static std::string to_string(const std::vector<string> &value)
    {
        return "{" + sbsearch::util::join(value, ",") + "}";
    }
};

#endif // SBSDB_POSTGRESQL_H_
