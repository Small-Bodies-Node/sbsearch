#ifndef SBS_EPHEMERIS_REMOVE_H_
#define SBS_EPHEMERIS_REMOVE_H_

#include "config.h"
#include "sbsearch.h"
#include "sbs-ephemeris/arguments.h"

using namespace sbsearch;

namespace sbsearch::sbs_ephemeris
{
    // remove the ephemeris points by target, optionally for a date range
    template <typename DB>
    void remove(const Arguments &args, SBSearch<DB> &sbs);
}

#endif