#ifndef SBS_EPHEMERIS_GET_H_
#define SBS_EPHEMERIS_GET_H_

#include <string>

#include "config.h"
#include "date.h"
#include "ephemeris.h"
#include "moving_target.h"
#include "sbsearch.h"
#include "sbs-ephemeris/arguments.h"

using namespace sbsearch;
using std::string;

namespace sbsearch::sbs_ephemeris
{
    // cli controller for sbs-ephemeris get action
    template <typename DB>
    void get(const Arguments &args, SBSearch<DB> &sbs);

    // Ephemeris segment length test for step_size = "auto".
    bool too_long(const Ephemeris::Datum &a, const Ephemeris::Datum &b);

    // get ephemeris data from horizons
    Ephemeris::Data get_from_horizons(const MovingTarget &target,
                                      string_view observer,
                                      const Date &start_date,
                                      const Date &stop_date,
                                      string_view step_size,
                                      const bool cache,
                                      const int recursion_step = 0,
                                      const bool refine = false);

    // Refine ephemeris data, replacing large steps with smaller steps.  Large
    // steps are identified with `too_long`.
    Ephemeris::Data refine_ephemeris(const MovingTarget &target,
                                     string_view observer,
                                     const Ephemeris::Data &data,
                                     const bool cache,
                                     const int recursion_step);
}

#endif