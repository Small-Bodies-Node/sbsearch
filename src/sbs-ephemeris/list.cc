#include <algorithm>
#include <iostream>
#include <utility>
#include <boost/json.hpp>

#include "config.h"
#include "json.h"
#include "moving_target.h"
#include "observatory.h"
#include "sbsearch.h"
#include "ephemeris/ephemeris.h"
#include "ephemeris/parallax_offset.h"
#include "sbsdb/get.h"
#include "sbsdb/postgresql.h"
#include "sbs-ephemeris/arguments.h"
#include "sbs-ephemeris/interpolate.h"

namespace json = boost::json;
using namespace sbsearch;

using sbsearch::ephemeris::Ephemeris;
using std::cerr;
using std::cout;
using std::endl;
using std::pair;

namespace sbsearch::sbs_ephemeris
{
    template <typename DB>
    void list(const Arguments &args, SBSearch<DB> &sbs)
    {
        MovingTarget target = sbsdb::get::moving_target(sbs.db(), args.target, !args.major_body);
        if (!target.moving_target_id())
            throw MovingTargetError("Target not found.");

        Ephemeris eph;
        if (args.interpolate > 0)
        {
            // For interpolation, get the whole ephemeris from the database and then
            // interpolate
            eph = sbsdb::get::ephemeris(sbs.db(), target);
            const Date start = args.start_date.value_or(Date(eph.data(0).mjd.value()));
            const Date stop = args.stop_date.value_or(Date(eph.data(-1).mjd.value())).mjd();
            eph = Ephemeris(eph.target(), interpolate(eph, start, stop, args.interpolate));
        }
        else
        {
            // If we are not interpolating, then search the database based on
            // start/stop dates
            eph = sbsdb::get::ephemeris(
                sbs.db(),
                target,
                args.start_date.value_or(0).mjd(),
                args.stop_date.value_or(100'000).mjd());
        }

        if (args.observer != "500@399")
        {
            Observatory observatory = sbsdb::get::observatory(sbs.db(), args.observer);
            eph = ephemeris::parallax_offset(eph, observatory);
        }

        std::ostream *os;
        std::ofstream outf;
        if (args.output_filename.empty())
            os = &cout;
        else
        {
            outf.open(args.output_filename, std::ios::trunc);
            os = &outf;
        }

        // output format
        eph.format.date = args.date_format;
        if (args.output_format == TABLE)
            *os << eph;
        else
            *os << json::value_from(eph);

        // close file stream
        if (outf.is_open())
            outf.close();
    }

    template void list(const Arguments &, SBSearch<sbsdb::Postgresql> &);
}