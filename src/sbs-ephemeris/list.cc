#include <algorithm>
#include <iostream>
#include <boost/json.hpp>

#include "config.h"
#include "ephemeris.h"
#include "moving_target.h"
#include "observatory.h"
#include "sbsearch.h"
#include "sbsdb/get.h"
#include "sbsdb/postgresql.h"
#include "sbs-ephemeris/arguments.h"

namespace json = boost::json;
using namespace sbsearch;
using std::cerr;
using std::cout;
using std::endl;

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

            // never extrapolate
            double start_mjd = std::max(args.start_date.value_or(Date(0)).mjd(), eph.data(0).mjd.value());
            double stop_mjd = std::min(args.stop_date.value_or(Date(100'000)).mjd(), eph.data(-1).mjd.value());

            Ephemeris interpolated;
            interpolated.target(target);
            for (double mjd = start_mjd; mjd < stop_mjd; mjd += args.interpolate)
                interpolated.append(eph.interpolate(mjd));

            interpolated.append(eph.interpolate(stop_mjd));
            eph = interpolated;
        }
        else
        {
            // If we are not interpolating, then search the database based on
            // start/stop dates
            cerr << "asdf\n";
            eph = sbsdb::get::ephemeris(
                sbs.db(),
                target,
                args.start_date.value_or(0).mjd(),
                args.stop_date.value_or(100'000).mjd());
            cerr << "asdf\n";
        }

        if (args.observer != "500@399")
        {
            Observatory observatory = sbsdb::get::observatory(sbs.db(), args.observer);
            eph = eph.parallax_offset(observatory);
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
        {
            json::array array = eph.as_json();
            *os << array;
        }

        // close file stream
        if (outf.is_open())
            outf.close();
    }

    template void list(const Arguments &, SBSearch<sbsdb::Postgresql> &);
}