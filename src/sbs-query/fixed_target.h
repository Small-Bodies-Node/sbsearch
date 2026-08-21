#ifndef SBS_SBSQUERY_FIXED_TARGET_H_
#define SBS_SBSQUERY_FIXED_TARGET_H_

#include <string>
#include "config.h"
#include "observation.h"
#include "sbsearch.h"
#include "sbs-query/arguments.h"

using namespace sbsearch;
using std::string;

namespace sbsearch::sbs_query::fixed_target
{
    template <typename DB>
    Observations query(const vector<string> &targets,
                       const Arguments &args,
                       SBSearch<DB> &sbs,
                       std::ostream *console);

    template <typename DB>
    Observations from_coordinates(string_view coordinates,
                                  const SearchOptions &search_options,
                                  SBSearch<DB> &sbs);
}

#endif
