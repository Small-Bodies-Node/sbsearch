#include <iostream>
#include <utility>

#include "config.h"
#include "date.h"
#include "ephemeris.h"
#include "sbs-ephemeris/interpolate.h"

using namespace sbsearch;
using std::cout;
using std::endl;
using std::pair;

namespace sbsearch::sbs_ephemeris
{
    Ephemeris::Data interpolate(Ephemeris &eph, const Date start, const Date stop, const double step)
    {
        // never extrapolate
        pair<double, double> eph_range = {
            eph.data(0).mjd.value(),
            eph.data(-1).mjd.value()};
        pair<double, double> args_range = {
            start.mjd(),
            stop.mjd()};
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
                 << interp_range.first << endl
                 << endl;

        if (args_range.second > eph_range.second)
            cout << "Interpolation beyond ephemeris not supported, stop date set to "
                 << interp_range.second << endl
                 << endl;

        Ephemeris::Data interpolated;
        for (double mjd = interp_range.first; mjd <= interp_range.second; mjd += step)
            interpolated.push_back(eph.interpolate(mjd));

        // always include the stop point
        if (std::fabs(interpolated.back().mjd.value() - stop.mjd()) > 1.0 / 86400.0)
            interpolated.push_back(eph.interpolate(stop.mjd()));

        return interpolated;
    }
}