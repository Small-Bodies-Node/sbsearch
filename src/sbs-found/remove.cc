#include <iostream>
#include <string>
#include <vector>

#include "config.h"
#include "found.h"
#include "moving_target.h"
#include "sbsdb.h"
#include "sbsearch.h"
#include "sbs-found/remove.h"

using namespace sbsearch;

using std::cout;
using std::endl;
using std::vector;

namespace sbsearch::sbs_found
{
    template <typename DB>
    void remove_found(SBSearch<DB> &sbs,
                      const double start_mjd,
                      const double stop_mjd)
    {
        int count = sbsdb::remove::found(sbs.db(), start_mjd, stop_mjd);
        cout << count << " found rows removed." << endl;
    }

    template <typename DB>
    void remove_found(SBSearch<DB> &sbs,
                      const vector<MovingTarget> &targets,
                      const double start_mjd,
                      const double stop_mjd)
    {
        for (auto const target : targets)
        {
            int count = sbsdb::remove::found(sbs.db(), target, start_mjd, stop_mjd);
            cout << count << " found rows removed." << endl;
        }
    }

    template void remove_found(SBSearch<sbsdb::Postgresql> &, const double, const double);
    template void remove_found(SBSearch<sbsdb::Postgresql> &, const vector<MovingTarget> &, const double, const double);
}