#ifndef SBS_EPHEMERIS_LIST_H_
#define SBS_EPHEMERIS_LIST_H_

#include "config.h"
#include "sbsearch.h"
#include "sbs-ephemeris/arguments.h"

using namespace sbsearch;

namespace sbsearch::sbs_ephemeris
{
    // list ephemeris data in the database, optionally for a date range,
    // optionally interpolating, optionally for an observatory
    template <typename DB>
    void list(const Arguments &args, SBSearch<DB> &sbs);
}

#endif