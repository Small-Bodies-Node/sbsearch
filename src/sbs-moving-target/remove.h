#ifndef SBS_MOVING_TARGET_REMOVE_H_
#define SBS_MOVING_TARGET_REMOVE_H_

#include <string_view>

#include "config.h"
#include "sbsearch.h"

using namespace sbsearch;
using std::string_view;

namespace sbsearch::sbs_moving_target
{
    template <typename DB>
    void remove(string_view name,
                const bool major_body,
                const bool force_remove,
                SBSearch<DB> &sbs);
}

#endif