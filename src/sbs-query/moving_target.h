#ifndef SBS_QUERY_MOVING_TARGET_H_
#define SBS_QUERY_MOVING_TARGET_H_

#include <string>
#include "config.h"
#include "found.h"
#include "observation.h"
#include "sbsearch.h"
#include "sbs-query/arguments.h"

using namespace sbsearch;
using std::string;

namespace sbsearch::sbs_query
{
    template <typename DB>
    const Founds query_moving_target(const Arguments &args, const string &designation, SBSearch<DB> &sbs);
}

#endif
