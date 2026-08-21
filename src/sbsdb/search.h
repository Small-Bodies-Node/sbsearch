#ifndef SBSDB_SEARCH_H_
#define SBSDB_SEARCH_H_

#include <cinttypes>
#include <optional>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "observation.h"

using std::optional;
using std::string;
using std::vector;

namespace sbsearch::sbsdb::search
{
    struct Options
    {
        double mjd_start = 0;
        double mjd_stop = 100000; // default: effectively search over all time
        optional<string> source;  // default: search all sources
    };

    // Search for observations matching the query terms, saving them to a
    // temporary database table, and returning the length of the table after
    // searching.
    //
    // The table is named search_observaitons_results and must be explicitly
    // dropped when no longer needed.  Get the results with results() and this
    // will be done automatically.
    //
    // In psql and sqlite3, it will only be accessible in from the given
    //  database connection.
    template <typename DB>
    int observations(
        DB *db, const vector<string> &query_terms, const Options &options = Options());

    // Get the results from observations() and remove the temporary table.
    template <typename DB>
    Observations results(DB *db);

    // Search for no more than this many terms at once
    const size_t maximum_query_terms = 1000;
}

#endif // SBSDB_SEARCH_H_
