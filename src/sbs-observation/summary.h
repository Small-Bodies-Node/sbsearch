#ifndef SBS_OBSERVATION_SUMMARY_H_
#define SBS_OBSERVATION_SUMMARY_H_

#include "./arguments.h"
#include "sbsearch.h"

namespace sbsearch::sbs_observation
{
    // generate a summary of observation coverage over the date range
    template <typename DB>
    void summary(const Arguments &args, SBSearch<DB> &sbs);
}

#endif