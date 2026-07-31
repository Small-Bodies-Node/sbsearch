#include <algorithm>
#include <iostream>
#include <utility>
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

            // never extrapolate
            pair<double, double> eph_range = {
                eph.data(0).mjd.value(),
                eph.data(-1).mjd.value()};
            pair<double, double> args_range = {
                args.start_date.value_or(Date(eph_range.first)).mjd(),
                args.stop_date.value_or(Date(eph_range.second)).mjd()};
            pair<double, double> interp_range = {
                std::max(args_range.first, eph_range.first),
                std::min(args_range.second, eph_range.second)};

            // range validation
            if (interp_range.first >= interp_range.second)
                throw std::runtime_error("Reqested time period is outside of the ephemeris range: " +
                                         std::to_string(eph_range.first) + " to " +
                                         std::to_string(eph_range.second));

            if (args_range.first < eph_range.first)
                cout << "Interpolation beyond ephemeris not supported, start date set to "
                     << interp_range.first << endl;

            if (args_range.second > eph_range.second)
                cout << "Interpolation beyond ephemeris not supported, stop date set to "
                     << interp_range.second << endl;

            Ephemeris interpolated;
            interpolated.target(target);
            for (double mjd = interp_range.first; mjd <= interp_range.second; mjd += args.interpolate)
                interpolated.append(eph.interpolate(mjd));

            eph = interpolated;
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