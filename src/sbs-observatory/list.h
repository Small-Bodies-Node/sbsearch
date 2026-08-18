#ifndef SBS_SBSOBSERVATORY_LIST_H_
#define SBS_SBSOBSERVATORY_LIST_H_

#include <string>
#include <string_view>

#include "config.h"
#include "sbsearch.h"

using namespace sbsearch;

using std::string;
using std::string_view;

namespace sbsearch::sbs_observatory
{
    template <typename DB>
    void list(string_view output_filename, SBSearch<DB> &sbs);
}

#endif