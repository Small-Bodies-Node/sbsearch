#include <iostream>
#include <iterator>
#include <string>
#include <s2/s1angle.h>

#include "config.h"
#include "ephemeris.h"
#include "files.h"
#include "horizons.h"
#include "moving_target.h"
#include "sbsearch.h"
#include "sbsdb/add.h"
#include "sbsdb/get.h"
#include "sbsdb/postgresql.h"
#include "sbs-ephemeris/get.h"

#define MAX_RECURSION 3

using namespace sbsearch;
using std::cerr;
using std::cout;
using std::endl;
using std::string;

namespace sbsearch::sbs_ephemeris
{
    template <typename DB>
    void get(const Arguments &args, SBSearch<DB> &sbs)
    {
        // default time range is from the database start/stop dates
        auto range = sbsdb::get::observations_date_range(sbs.db());
        if ((args.action != "remove") && (!args.start_date || !args.stop_date) && (!range.first || !range.second))
            throw EphemerisError("Observations database is empty: --start and --stop are required.");

        const Date start_date = args.start_date ? args.start_date.value() : range.first.value();
        const Date stop_date = args.stop_date ? args.stop_date.value() : range.second.value();

        MovingTarget target(args.target, !args.major_body);

        cout << "Fetching ephemeris from Horizons for " << target.designation() << " from "
             << start_date.iso() << " to " << stop_date.iso()
             << " with step size " << args.time_step
             << "." << endl;

        Ephemeris eph(target, get_from_horizons(target, args.observer, start_date, stop_date, args.time_step, args.cache));
        cout << endl;

        // write to screen or file?
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

    bool too_long(const Ephemeris::Datum &a, const Ephemeris::Datum &b)
    {
        S1Angle angle(a.as_s2point(), b.as_s2point());
        return angle.degrees() > 0.5;
    }

    Ephemeris::Data get_from_horizons(const MovingTarget &target,
                                      const string &observer,
                                      const Date &start_date,
                                      const Date &stop_date,
                                      const string &time_step,
                                      const bool cache,
                                      const int recursion_step,
                                      const bool refine)
    {
        // fail safe, code below should never reach here
        if (recursion_step > MAX_RECURSION)
            throw(std::runtime_error("Maximum recursion on get_from_horizons exceeded."));

        // setup time periods to request; on the first run, split up
        // start_date/stop_date, but on successive runs just use them as is
        vector<std::pair<Date, Date>> dates;
        if (recursion_step == 0)
        {
            // request a maximum of 1 year at a time for VAR steps, 5 years otherwise
            const bool variable_steps = (time_step.find("VAR") == string::npos) || (time_step == "auto");
            const int years = variable_steps ? 1 : 5;
            dates = date_ranges(start_date, stop_date, years * 365.);
        }
        else
            dates.emplace_back(start_date, stop_date);

        // iterate over time periods
        bool first = true;
        Ephemeris::Data data;
        for (auto const &[start, stop] : dates)
        {
            cout << string(recursion_step + 1, '>')
                 << " Current step " << start.iso()
                 << " to " << stop.iso()
                 << " with step " << time_step
                 << "." << endl;

            Horizons horizons(target,
                              observer,
                              start,
                              stop,
                              time_step == "auto" ? "VAR 900" : time_step,
                              cache);
            Ephemeris::Data new_data = horizons.get_ephemeris_data();

            // verify the last step when using VAR steps
            if ((time_step == "auto" || time_step.find("VAR") != string::npos) &&
                (std::fabs(new_data.back().mjd.value() - stop.mjd()) > 0.0001))
            {
                Logger::warning() << "Variable time step detected, but end date was not returned.  Retrying." << endl;
                horizons.time_step("1");
                new_data.push_back(horizons.get_ephemeris_data().back());
            }

            // identify large steps and replace with finer steps
            if ((refine || time_step == "auto") && recursion_step < MAX_RECURSION)
                new_data = refine_ephemeris(target, observer, new_data, cache, recursion_step + 1);

            // skip the first point if it was already added on the last step
            auto copy_from = new_data.cbegin();
            if (!first)
                std::advance(copy_from, 1);

            std::copy(copy_from, new_data.cend(), std::back_inserter(data));
            first = false;
        }

        return data;
    }

    Ephemeris::Data refine_ephemeris(const MovingTarget &target,
                                     const string &observer,
                                     const Ephemeris::Data &data,
                                     const bool cache,
                                     const int recursion_step)
    {
        // results
        Ephemeris::Data refined;

        // initialize iterators
        auto refine_start = data.cbegin();
        auto refine_stop = data.cbegin();

        // loop over all the data looking for segments to refine
        int n;
        Date refine_start_date, refine_stop_date;
        do
        {
            // identify the next step to refine
            for (refine_start = refine_stop; refine_start < std::prev(data.cend()); refine_start++)
                if (too_long(*refine_start, *std::next(refine_start)))
                    break;

            // save the good data since the last refinement, up to the last point
            std::copy(refine_stop, refine_start, std::back_inserter(refined));

            // if we are at the end of the new data, add the last point and we're done
            if (refine_start >= std::prev(data.cend()))
            {
                refined.push_back(*refine_start);
                break;
            }

            // identify the end of the step to refine
            for (refine_stop = std::next(refine_start); refine_stop < std::prev(data.cend()); std::advance(refine_stop, 1))
                if (!too_long(*refine_stop, *std::next(refine_stop)))
                    break;

            // target 1/4 degree per step
            double length = 0;
            for (auto it = refine_start; it < refine_stop; it++)
                length += S1Angle((*it).as_s2point(), (*std::next(it)).as_s2point()).degrees();

            n = 4 * static_cast<int>(std::ceil(length));
            Logger::debug() << "Refining " << length << " deg arc with " << n << " steps." << endl;
            refine_start_date = (*refine_start).mjd.value();
            refine_stop_date = (*refine_stop).mjd.value();

            // get and save the data
            Ephemeris::Data refined_data = get_from_horizons(target,
                                                             observer,
                                                             refine_start_date,
                                                             refine_stop_date,
                                                             std::to_string(n),
                                                             cache,
                                                             recursion_step,
                                                             true);

            // again, append up to the last point
            std::copy(refined_data.cbegin(),
                      std::prev(refined_data.cend()),
                      std::back_inserter(refined));

            // if we are at the end, append the last point
            if (refine_stop >= std::prev(data.cend()))
                refined.push_back(refined_data.back());
        } while (refine_stop < std::prev(data.cend()));

        return refined;
    }

    template void get(const Arguments &, SBSearch<sbsdb::Postgresql> &);
}