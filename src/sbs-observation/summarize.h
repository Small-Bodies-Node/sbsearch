#ifndef SBS_SBSOBSERVATION_SUMMARY_H_
#define SBS_SBSOBSERVATION_SUMMARY_H_

#include "sbs-observation/arguments.h"
#include "sbsearch.h"

namespace sbsearch::sbs_observation
{
    // Summarize observation coverage over the date range
    template <typename DB>
    void summarize(const Arguments &args, SBSearch<DB> &sbs);
}

#endif