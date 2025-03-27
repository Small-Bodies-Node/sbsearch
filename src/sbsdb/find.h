#ifndef SBSDB_FIND_H_
#define SBSDB_FIND_H_

#include <cinttypes>
#include <optional>
#include <string>
#include <unordered_set>
#include <vector>

using std::optional;
using std::string;
using std::vector;

namespace sbsearch::sbsdb::find
{
    struct Options
    {
        double mjd_start = 0;
        double mjd_stop = 100000; // default: effectively search over all time
        optional<string> source;  // default: search all sources
    };

    template <typename DB>
    std::unordered_set<int64_t> observations(DB &db,
                                             const vector<string> &query_terms,
                                             const Options &options = Options());

    // Search for no more than this many terms at once
    const size_t maximum_query_terms = 100;
}

#endif // SBSDB_FIND_H_
