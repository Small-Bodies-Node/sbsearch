#ifndef SBSDB_POSTGRESQL_H_
#define SBSDB_POSTGRESQL_H_

#include <optional>
#include <string>
#include <pqxx/pqxx>

#include "sbsdb.h"
#include "../ephemeris.h"
#include "../exceptions.h"
#include "../found.h"
#include "../moving_target.h"
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
         * @return std::vector<T>
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
                                std::vector<T>>
        {
            const int nargs = sizeof...(args);

            pqxx::result result;

            if (nargs == 0)
                result = work_.exec(statement);
            else
                result = work_.exec_params(statement, pqxx::params(args...));

            std::vector<T> v(result.size());
            std::transform(result.begin(), result.end(), v.begin(),
                           [&](const pqxx::row &row)
                           { return row_as<T>(row); });

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
        T row_as(const pqxx::row &row);
    };
}

#endif // SBSDB_POSTGRESQL_H_
