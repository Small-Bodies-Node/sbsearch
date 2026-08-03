#ifndef SBS_FOUND_LIST_H_
#define SBS_FOUND_LIST_H_

#include <string>
#include <vector>

#include "config.h"
#include "cli.h"
#include "found.h"
#include "moving_target.h"
#include "sbsearch.h"

using std::string;
using std::vector;

namespace sbsearch::sbs_found
{
    // Remove found results not in source list.  Does nothing if source list is
    // empty.
    void filter_sources(const vector<string> &sources, Founds &founds);

    // List observations of any moving target
    template <typename DB>
    void list_found(SBSearch<DB> &sbs,
                    const double start_mjd,
                    const double stop_mjd,
                    const vector<string> &sources,
                    const string output_filename,
                    const cli::OutputFormat output_format);

    // List observations of specific moving targets
    template <typename DB>
    void list_found(SBSearch<DB> &sbs,
                    const vector<MovingTarget> &targets,
                    const double start_mjd,
                    const double stop_mjd,
                    const vector<string> &sources,
                    const string output_filename,
                    const cli::OutputFormat output_format);

    // List a set of found results.
    void list_found(const Founds &founds,
                    const string output_filename,
                    const cli::OutputFormat output_format);
}
#endif
