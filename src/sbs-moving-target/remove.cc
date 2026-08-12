#include <iostream>
#include <string>
#include <vector>

#include "cli.h"
#include "config.h"
#include "date.h"
#include "moving_target.h"
#include "sbs-moving-target/remove.h"
#include "sbsdb.h"
#include "sbsearch.h"

using namespace sbsearch;
using sbsearch::SBSearch;
using sbsearch::cli::confirm;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace sbsearch::sbs_moving_target
{
    template <typename DB>
    void remove(string_view name,
                const bool major_body,
                const bool force_remove,
                SBSearch<DB> &sbs)
    {
        MovingTarget target = sbsdb::get::moving_target(sbs.db(), name, !major_body);
        if (!target.moving_target_id())
            cout << name << " not in the database.\n";
        else
        {
            if (force_remove || confirm("Remove target " + target.to_string() + ", ephemeris, and found observations?"))
            {
                int count = sbsdb::remove::found(sbs.db(), target, 0, 100'000);
                cout << "Removed " << count << " found entries." << endl;

                count = sbsdb::remove::ephemeris(sbs.db(), target);
                cout << "Removed " << count << " ephemeris points." << endl;

                sbsdb::remove::moving_target(sbs.db(), target);
                cout << "Removed target " << target << "." << endl;
            }
        }
    }

    template void remove(string_view, const bool, const bool, SBSearch<sbsdb::Postgresql> &);
}
