#ifndef SBS_MOVING_TARGET_ADD_H_
#define SBS_MOVING_TARGET_ADD_H_

#include <string>
#include <string_view>
#include <vector>

#include "config.h"
#include "sbsearch.h"

using namespace sbsearch;

using std::string;
using std::string_view;
using std::vector;

namespace sbsearch::sbs_moving_target
{
    template <typename DB>
    void add(string_view name,
             const bool major_body,
             const vector<string> &alternate_names,
             SBSearch<DB> &sbs);
}

#endif