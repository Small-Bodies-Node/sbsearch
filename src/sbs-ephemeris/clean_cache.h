#ifndef SBS_SBSEPHEMERIS_CLEAN_CACHE_H_
#define SBS_SBSEPHEMERIS_CLEAN_CACHE_H_

#include "config.h"
#include "ephemeris/ephemeris.h"
#include "sbs-ephemeris/arguments.h"

using namespace sbsearch;
using sbsearch::ephemeris::Ephemeris;
using std::string;

namespace sbsearch::sbs_ephemeris
{
    // Remove files from the ephemeris data cache older than max_age days.
    void clean_cache(double max_age);
}

#endif