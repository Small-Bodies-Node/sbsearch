#include <algorithm>
#include <cinttypes>
#include <optional>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

#include "find.h"
#include "postgresql.h"
#include "../intersection.h"
#include "../observation.h"

using std::endl;
using std::string;
using std::vector;

namespace sbsearch::sbsdb::find
{
    template <typename DB>
    Observations observations(DB *db, const vector<string> &query_terms, const Options &options)
    {
        Logger::debug() << query_terms.size() << " query terms to search." << std::endl;

        Observations matches;
        const bool use_transaction = db->template begin();
        try
        {
            // Create a temporary table for the results.  The name is constant,
            // so we do need to explicitly drop it when finished.  It should not
            // be accessible to other connections (psql and sqlite behaviors).
            db->template execute("CREATE TEMPORARY TABLE find_observations_results (LIKE observations)");

            // Query database with terms, but not too many at once
            vector<string> subset;
            subset.reserve(maximum_query_terms);

            string statement = "INSERT INTO find_observations_results SELECT * FROM observations WHERE";
            if constexpr (std::is_same_v<DB, Postgresql> == true)
                statement += " terms && $1";
            else // Sqlite
                statement = " terms MATCH $1";

            if (options.source)
                statement += " AND source = $2 AND mjd_start >= $3 AND mjd_stop <= $4";
            else
                statement += " AND mjd_start >= $2 AND mjd_stop <= $3";

            // Store results in a temporary table
            for (int i = 0; i < query_terms.size(); i += maximum_query_terms)
            {
                const int j = std::min(query_terms.size(), i + maximum_query_terms);
                subset.assign(query_terms.begin() + i, query_terms.begin() + j);
                if (options.source)
                    db->template execute(statement, subset, options.source, options.mjd_start, options.mjd_stop);
                else
                    db->template execute(statement, subset, options.mjd_start, options.mjd_stop);
            }

            // Get the results
            matches = Observations(db->template get_all_observations("find_observations_results"));

            // Done with the temporary table
            db->template execute("DROP TABLE find_observations_results");
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

        matches.remove_duplicate_observation_ids();

        // Logger::info() << matches.size() << " approximate matches." << endl;

        return matches;
    };

    template Observations observations(Postgresql *, const vector<string> &, const Options &);
}