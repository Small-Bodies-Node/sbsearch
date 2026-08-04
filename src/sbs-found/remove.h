#ifndef SBS_FOUND_REMOVE_H_
#define SBS_FOUND_REMOVE_H_

#include <string>
#include <vector>

#include "config.h"
#include "found.h"
#include "moving_target.h"
#include "sbsearch.h"

using std::string;
using std::vector;

namespace sbsearch::sbs_found
{
    // Remove observations of all moving targets within date range
    template <typename DB>
    void remove_found(SBSearch<DB> &sbs,
                      const double start_mjd,
                      const double stop_mjd);

    // List observations of any moving target
    template <typename DB>
    void remove_found(SBSearch<DB> &sbs,
                      const vector<MovingTarget> &targets,
                      const double start_mjd,
                      const double stop_mjd);
}

#endif
