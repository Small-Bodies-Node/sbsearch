#include <algorithm>
#include <set>
#include <string>
#include <vector>

#include "config.h"
#include "found.h"

using namespace sbsearch;

using std::string;
using std::vector;

namespace sbsearch::sbs_found
{
    void filter_sources(const vector<string> &sources, Founds &founds)
    {
        if (sources.size() > 0)
        {
            std::set<string> source_set(sources.begin(), sources.end());

            auto is_requested_source = [&source_set](const Found &found)
            {
                return source_set.count(found.observation.source()) != 0;
            };

            founds.data.erase(std::remove_if(founds.data.begin(),
                                             founds.data.end(),
                                             is_requested_source),
                              founds.data.end());
        }
    }
}