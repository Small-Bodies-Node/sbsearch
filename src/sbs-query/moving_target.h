#ifndef SBS_SBSQUERY_MOVING_TARGET_H_
#define SBS_SBSQUERY_MOVING_TARGET_H_

#include <string>
#include "config.h"
#include "found.h"
#include "observation.h"
#include "sbsearch.h"
#include "sbs-query/arguments.h"

using namespace sbsearch;
using std::string;

namespace sbsearch::sbs_query::moving_target
{
    template <typename DB>
    Founds query(const vector<string> &target_names,
                 const Arguments &args,
                 SBSearch<DB> &sbs,
                 std::ostream *console);

    template <typename DB>
    Founds from_ephemeris(Ephemeris &eph,
                          const vector<string> &sources,
                          SearchOptions &search_options,
                          SBSearch<DB> &sbs,
                          std::ostream *console);

    template <typename DB>
    Founds from_database(string_view name,
                         const Date &eph_start_date,
                         const Date &eph_stop_date,
                         const vector<string> &sources,
                         SearchOptions &search_options,
                         SBSearch<DB> &sbs,
                         std::ostream *console);

    template <typename DB>
    Founds from_ephemeris_file(string_view name,
                               string_view eph_file,
                               const vector<string> &sources,
                               SearchOptions &search_options,
                               SBSearch<DB> &sbs,
                               std::ostream *console);

    template <typename DB>
    Founds from_orbit_file(string_view name,
                           string_view orbit_file,
                           const Date &eph_start_date,
                           const Date &eph_stop_date,
                           string_view step_size,
                           bool cache,
                           const vector<string> &sources,
                           SearchOptions &search_options,
                           SBSearch<DB> &sbs,
                           std::ostream *console);

    template <typename DB>
    Founds from_horizons(const MovingTarget &target,
                         const Date &eph_start_date,
                         const Date &eph_stop_date,
                         string_view step_size,
                         const bool cache,
                         const vector<string> &sources,
                         SearchOptions &search_options,
                         SBSearch<DB> &sbs,
                         std::ostream *console);
}

#endif
