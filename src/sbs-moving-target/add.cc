#include <iostream>
#include <string>
#include <vector>

#include "config.h"
#include "logging.h"
#include "moving_target.h"
#include "sbs-moving-target/add.h"
#include "sbsdb.h"
#include "sbsearch.h"

using sbsearch::SBSearch;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::string_view;
using std::vector;

using namespace sbsearch;

namespace sbsearch::sbs_moving_target
{
    template <typename DB>
    void add(string_view name,
             const bool major_body,
             const vector<string> &alternate_names,
             SBSearch<DB> &sbs)
    {
        MovingTarget target{name, !major_body};
        target.add_names(alternate_names.begin(), alternate_names.end());
        sbsdb::add::moving_target(sbs.db(), target);
        cout << "Added " << target << "\n";
    }

    template void add(string_view, const bool, const vector<string> &, SBSearch<sbsdb::Postgresql> &);
}
