#ifndef SBS_EPHEMERIS_LIST_H_
#define SBS_EPHEMERIS_LIST_H_

#include "config.h"
#include "sbsearch.h"
#include "sbs-ephemeris/arguments.h"

using namespace sbsearch;

namespace sbsearch::sbs_ephemeris
{
    // cli controller for sbs-ephemeris list action
    template <typename DB>
    void list(const Arguments &args, SBSearch<DB> &sbs);
}

#endif