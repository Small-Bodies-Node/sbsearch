#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <boost/program_options.hpp>
#include <curl/curl.h>

#include "config.h"
#include "exceptions.h"
#include "moving_target.h"
#include "sbsearch.h"
#include "sbsdb/get.h"
#include "sbsdb/postgresql.h"
#include "sbsdb/remove.h"
#include "sbs-ephemeris/arguments.h"

using namespace sbsearch;

using sbsearch::ephemeris::Ephemeris;
using std::cerr;
using std::cout;
using std::endl;
using std::string;

namespace sbsearch::sbs_ephemeris
{
    template <typename DB>
    void remove(const Arguments &args, SBSearch<DB> &sbs)
    {
        MovingTarget target = sbsdb::get::moving_target(sbs.db(), args.target, !args.major_body);
        if (!target.moving_target_id())
            throw MovingTargetError(args.target + " not in database");

        int count;
        if (args.remove_all)
            count = sbsdb::remove::ephemeris(sbs.db(), target);
        else
            count = sbsdb::remove::ephemeris(
                sbs.db(),
                target,
                args.start_date.value_or(Date(0)).mjd(),
                args.stop_date.value_or(Date(100'000)).mjd());

        cout << count << " ephemeris rows removed." << endl;
    }

    template void remove(const Arguments &, SBSearch<sbsdb::Postgresql> &);
}