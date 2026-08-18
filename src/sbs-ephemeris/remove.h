#ifndef SBS_SBSEPHEMERIS_REMOVE_H_
#define SBS_SBSEPHEMERIS_REMOVE_H_

#include "config.h"
#include "sbsearch.h"
#include "sbs-ephemeris/arguments.h"

using namespace sbsearch;

namespace sbsearch::sbs_ephemeris
{
    // cli controller for sbs-ephemeris remove action
    template <typename DB>
    void remove(const Arguments &args, SBSearch<DB> &sbs);
}

#endif