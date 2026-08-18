#ifndef SBS_SBSFOUND_FILTER_SOURCES_H_
#define SBS_SBSFOUND_FILTER_SOURCES_H_

#include <string>
#include <vector>

#include "config.h"
#include "found.h"

using std::string;
using std::vector;

namespace sbsearch::sbs_found
{
    // Remove found results not in source list.  Does nothing if source list is
    // empty.
    void filter_sources(const vector<string> &sources, Founds &founds);
}

#endif
