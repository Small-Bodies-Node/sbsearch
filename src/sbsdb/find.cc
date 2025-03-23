#include <algorithm>
#include <cinttypes>
#include <optional>
#include <type_traits>
#include <unordered_set>
#include <vector>

#include "find.h"
#include "postgresql.h"
// #include "../ephemeris.h"
// #include "../exceptions.h"
// #include "../found.h"
// #include "../moving_target.h"
// #include "../observation.h"
// #include "../observatory.h"

using std::endl;
using std::string;
using std::vector;

namespace sbsearch::sbsdb::find
{
    template <typename DB>
    std::unordered_set<int64_t> observations(DB &db, const vector<string> &query_terms, const Options &options)
    {
        Logger::debug() << query_terms.size() << " query terms to search." << std::endl;

        // Query database with terms, but not too many at once
        vector<string> subset;
        subset.reserve(maximum_query_terms);

        string statement = "SELECT observation_id FROM observations WHERE";
        if constexpr (std::is_same_v<DB, Postgresql> == true)
            statement += " terms && $1";
        else // Sqlite
            statement = " terms MATCH $1";

        int index = 1;
        if (options.source)
            statement += " AND source = $2 AND mjd_start >= $3 AND mjd_stop <= $4";
        else
            statement += " AND mjd_start >= $2 AND mjd_stop <= $3";

        std::unordered_set<int64_t> approximate_matches;
        for (int i = 0; i < query_terms.size(); i += maximum_query_terms)
        {
            const int j = std::min(query_terms.size(), i + maximum_query_terms);
            subset.assign(query_terms.begin() + i, query_terms.begin() + j);

            // std::copy(query_terms.begin() + i,
            //           query_terms.begin() + j,
            //           std::ostream_iterator<string>(cerr, " "));

            vector<int64_t> result;
            if (options.source)
                result = db.template get_many<int64_t>(statement,
                                                       subset,
                                                       options.source,
                                                       options.mjd_start,
                                                       options.mjd_stop);
            else
                result = db.template get_many<int64_t>(statement,
                                                       subset,
                                                       options.mjd_start,
                                                       options.mjd_stop);

            for (auto const &observation_id : result)
                approximate_matches.insert(observation_id);

            Logger::debug() << "Searched " << j << " of "
                            << query_terms.size() << " query terms, found "
                            << result.size() << " approximate matches."
                            << endl;
        }

        return approximate_matches;
    };

    template std::unordered_set<int64_t> observations(Postgresql &, const vector<string> &, const Options &);
}