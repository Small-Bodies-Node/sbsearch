#ifndef SBS_SBSEPHEMERIS_ADD_H_
#define SBS_SBSEPHEMERIS_ADD_H_

#include "config.h"
#include "date.h"
#include "moving_target.h"
#include "sbsearch.h"
#include "ephemeris/ephemeris.h"
#include "sbs-ephemeris/arguments.h"

using namespace sbsearch;
using sbsearch::ephemeris::Ephemeris;
using std::string;

namespace sbsearch::sbs_ephemeris
{
    // cli controller for sbs-ephemeris add action
    template <typename DB>
    void add(const Arguments &args, SBSearch<DB> &sbs);

    // add ephemeris data from file
    template <typename DB>
    int add_from_file(string_view file, MovingTarget &target, SBSearch<DB> &sbs);

    // add ephemeris data from horizons
    template <typename DB>
    int add_from_horizons(const MovingTarget &target,
                          const Date &start_date,
                          const Date &stop_date,
                          string_view step_size,
                          bool cache,
                          SBSearch<DB> &sbs);
}

#endif