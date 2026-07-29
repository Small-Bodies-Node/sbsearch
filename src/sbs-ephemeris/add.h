#ifndef SBS_EPHEMERIS_ADD_H_
#define SBS_EPHEMERIS_ADD_H_

#include "config.h"
#include "sbsearch.h"
#include "sbs-ephemeris/arguments.h"

using namespace sbsearch;

namespace sbsearch::sbs_ephemeris
{
    // add ephemeris data from file or horizons
    template <typename DB>
    void add(const Arguments &args, SBSearch<DB> &sbs);
}

#endif