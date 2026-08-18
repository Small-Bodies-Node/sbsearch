#include <iostream>
#include <string>
#include <vector>

#include "config.h"
#include "moving_target.h"
#include "sbsdb.h"
#include "sbsearch.h"
#include "table.h"
#include "sbs-moving-target/list.h"
#include "util/string.h"

using namespace sbsearch;
using namespace sbsearch::table;
using sbsearch::SBSearch;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace sbsearch::sbs_moving_target
{
    template <typename DB>
    void list(SBSearch<DB> &sbs)
    {
        using sbsearch::util::join;

        vector<string> designations, alternates;
        for (const MovingTarget &target : sbsdb::get::all_moving_targets(sbs.db()))
        {
            designations.emplace_back(target.designation());
            auto alt = target.alternate_names();
            alternates.push_back(join<string>({alt.begin(), alt.end()}, ", "));
        }

        table::Table tab;
        tab.add(Column("designation", "%s", designations));
        tab.add(Column("alternates", "%s", alternates));
        std::cout << tab;
    }

    template void list(SBSearch<sbsdb::Postgresql> &);
}
