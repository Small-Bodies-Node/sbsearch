#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <iostream>
#include <optional>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

#include "config.h"
#include "intersection.h"
#include "logging.h"
#include "observation.h"
#include "sbsdb/search.h"
#include "sbsdb/postgresql.h"
#include "util/string.h"

using namespace std::chrono_literals;
using std::endl;
using std::string;
using std::vector;

namespace sbsearch::sbsdb::search
{
    template <typename DB>
    int observations(DB *db, const vector<string> &query_terms, const Options &options)
    {
        const bool use_transaction = db->template begin();
        try
        {
            // Create a temporary table for the results.  It should not be
            // accessible to other connections (psql and sqlite behaviors). We
            // are using "IF NOT EXISTS" so that multiple calls can append new
            // results.
            db->template execute("CREATE TEMPORARY TABLE IF NOT EXISTS search_observaitons_results (LIKE observations)");

            // Query database with terms, but not too many at once
            vector<string> query_subset;
            query_subset.reserve(maximum_query_terms);

            string statement = "INSERT INTO search_observaitons_results SELECT * FROM observations WHERE";
            if constexpr (std::is_same_v<DB, Postgresql> == true)
                statement += " terms && $1 AND mjd_range && numrange($2, $3)";
            else // Sqlite
                statement = " terms MATCH $1 AND mjd_start < $3 AND mjd_stop > $2";

            if (options.source)
                statement += " AND source = $4";

            // Store results in a temporary table
            auto start = std::chrono::steady_clock::now();
            for (int i = 0; i < query_terms.size(); i += maximum_query_terms)
            {
                const int j = std::min(query_terms.size(), i + maximum_query_terms);
                query_subset.assign(query_terms.begin() + i, query_terms.begin() + j);
                if (options.source)
                    db->template execute(statement,
                                         query_subset,
                                         options.mjd_start,
                                         options.mjd_stop,
                                         options.source);
                else
                    db->template execute(statement,
                                         query_subset,
                                         options.mjd_start,
                                         options.mjd_stop);

                if (options.debug)
                {
                    auto end = std::chrono::steady_clock::now();
                    const std::chrono::duration<double> diff = end - start;
                    if (diff > 300ms)
                        Logger::debug() << "Long query time detected (" << diff.count() << "s): "
                                        << statement << " "
                                        << "array['" << util::join(query_subset, "','") << "'] "
                                        << options.mjd_start << " "
                                        << options.mjd_stop << " "
                                        << options.source.value_or("all sources") << " "
                                        << std::endl;
                    start = end;
                }
            }
        }
        catch (const SBSException &err)
        {
            if (use_transaction)
                db->template rollback();
            throw err;
        }
        catch (const std::exception &err)
        {
            if (use_transaction)
                db->template rollback();
            Logger::error() << err.what() << endl;
            throw err;
        }

        int64_t n = db->template get_one<int>("SELECT COUNT(DISTINCT(observation_id)) FROM search_observaitons_results");
        Logger::debug() << "Searched for " << query_terms.size() << " query terms and collected " << n << " approximate matches." << endl;

        if (use_transaction)
            db->template commit();

        return n;
    };

    template <typename DB>
    Observations results(DB *db)
    {
        auto matches = Observations(db->template get_all_observations("search_observaitons_results"));

        // Done with the temporary table
        db->template execute("DROP TABLE search_observaitons_results");

        return matches;
    }

    template int observations(Postgresql *, const vector<string> &, const Options &);
    template Observations results(Postgresql *);
}