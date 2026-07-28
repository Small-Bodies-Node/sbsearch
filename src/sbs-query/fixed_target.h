#ifndef SBS_QUERY_FIXED_TARGET_H_
#define SBS_QUERY_FIXED_TARGET_H_

#include <string>
#include "config.h"
#include "observation.h"
#include "sbsearch.h"
#include "sbs-query/arguments.h"

using namespace sbsearch;
using std::string;

namespace sbsearch::sbs_query
{
    template <typename DB>
    const Observations query_fixed_target(const Arguments &args, const string &coordinates, SBSearch<DB> &sbs);
}

#endif
