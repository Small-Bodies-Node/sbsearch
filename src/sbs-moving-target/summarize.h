#ifndef SBS_SBSMOVING_TARGET_SUMMARIZE_H_
#define SBS_SBSMOVING_TARGET_SUMMARIZE_H_

#include <optional>

#include "config.h"
#include "date.h"
#include "sbsearch.h"

using namespace sbsearch;

using std::optional;

namespace sbsearch::sbs_moving_target
{
    template <typename DB>
    void summarize(const optional<Date> &start_date,
                   const optional<Date> &stop_date,
                   SBSearch<DB> &sbs);
}

#endif