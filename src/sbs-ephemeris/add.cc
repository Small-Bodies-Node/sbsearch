#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

#include "config.h"
#include "ephemeris.h"
#include "files.h"
#include "horizons.h"
#include "moving_target.h"
#include "sbsearch.h"
#include "sbsdb/add.h"
#include "sbsdb/get.h"
#include "sbsdb/postgresql.h"
#include "sbs-ephemeris/add.h"

using namespace sbsearch;
using std::cerr;
using std::cout;
using std::endl;
using std::string;

namespace sbsearch::sbs_ephemeris
{
    Ephemeris ephemeris_helper_(const MovingTarget &target,
                                const string &observer,
                                const Date &start,
                                const Date &stop,
                                const string &time_step,
                                const bool cache)
    {
        // initial step size for "auto" mode is 1 day, otherwise use the
        // requested step
        string step = time_step == "auto" ? "1 day" : time_step;

        Horizons horizons(target, observer, start, stop, step, cache);
        Ephemeris eph(target, horizons.get_ephemeris_data());

        if ((time_step.find("VAR") != string::npos) && (eph.data(-1).mjd < stop.mjd()))
        {
            // always verify the last item when using VAR steps
            Logger::warning() << "Variable time step detected, but end date was not returned.  Retrying." << endl;
            horizons.time_step("1");
            auto data = horizons.get_ephemeris_data();
            eph.append(Ephemeris(target, data)[1]);
        }
        else if (time_step == "auto")
        {
            // For "auto," recursively identify large steps and split them up
            Ephemeris::Data revised;
            bool first = true;
            for (const auto segment : eph.segments())
            {
                Ephemeris::Data new_data;
                const double length = segment.as_polyline().GetLength().degrees();
                if (length < 0.5)
                {
                    new_data = segment.data();
                }
                else
                {
                    step = std::to_string((static_cast<int>(std::ceil(4 * length))));
                    new_data = ephemeris_helper_(target,
                                                 observer,
                                                 segment.data(0).mjd.value(),
                                                 segment.data(1).mjd.value(),
                                                 step,
                                                 cache)
                                   .data();
                }

                auto begin = new_data.begin();
                if (!first)
                    // skip the first point (already added)
                    std::advance(begin, 1);

                std::copy(begin, new_data.end(), std::back_inserter(revised));
                first = false;
            }

            eph = Ephemeris(target, revised);
        }
        return eph;
    }

    template <typename DB>
    void add(const Arguments &args, SBSearch<DB> &sbs)
    {
        MovingTarget target = sbsdb::get::moving_target(sbs.db(), args.target, !args.major_body);
        std::cerr << target << std::endl;
        if (!target.moving_target_id())
        {
            sbsdb::add::moving_target(sbs.db(), target);
            cout << "Added moving target " << target.designation()
                 << " to the database with ID " << target.moving_target_id().value()
                 << "." << endl;
        }

        int count = 0;
        if (!args.file.empty())
        {
            cout << "Reading ephemeris from file " << args.file << ".\n";
            string table = read_file(args.file);
            Ephemeris eph = Ephemeris(target, Horizons::parse(table));
            sbs.add_ephemeris(eph);
            count = eph.num_vertices();
        }
        else
        {
            cout << "Fetching ephemeris from Horizons API for " << target.designation() << " from "
                 << args.start_date.value().iso() << " to " << args.stop_date.value().iso()
                 << "." << endl;

            // request a maximum of 1 year at a time for VAR steps, 5 years otherwise
            const int years = (args.time_step.find("VAR") == string::npos) ? 5 : 1;
            auto dates = date_ranges(args.start_date.value(), args.stop_date.value(), years * 365.);

            bool first = true;
            Ephemeris ephemeris(target, {});
            for (auto const &[start, stop] : dates)
            {
                Ephemeris eph = ephemeris_helper_(target,
                                                  args.observer,
                                                  start,
                                                  stop,
                                                  args.time_step,
                                                  args.cache);

                // trim the first point... it was already added on the last step
                if (!first)
                    ephemeris.append({std::next(eph.data().begin()), eph.data().end()});
                else
                    ephemeris.append(eph);

                first = false;
            }
            count = ephemeris.num_vertices();
            sbs.add_ephemeris(ephemeris);
        }

        if (count == 0)
            throw std::runtime_error("Empty ephemeris data.");
        else
            cout << "Added " << count << " ephemeris epochs.\n\n";
    }

    template void add(const Arguments &, SBSearch<sbsdb::Postgresql> &);
}