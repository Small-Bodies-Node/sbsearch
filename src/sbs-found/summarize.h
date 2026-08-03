#ifndef SBS_FOUND_SUMMARY_H_
#define SBS_FOUND_SUMMARY_H_

#include "config.h"
#include "cli.h"
#include "moving_target.h"
#include "sbsearch.h"

using namespace sbsearch;

using std::vector;

namespace sbsearch::sbs_found
{
    // Summarize observations of moving targets
    template <typename DB>
    void summarize_found(SBSearch<DB> &sbs,
                         const vector<MovingTarget> &targets,
                         const double start_mjd,
                         const double stop_mjd,
                         const vector<string> &sources,
                         const string output_filename,
                         const cli::OutputFormat output_format);
}

#endif
