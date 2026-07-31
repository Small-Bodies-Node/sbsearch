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
#include "sbs-ephemeris/arguments.h"
#include "sbs-ephemeris/get.h"

using namespace sbsearch;
using std::cerr;
using std::cout;
using std::endl;
using std::string;

namespace sbsearch::sbs_ephemeris
{
    template <typename DB>
    void add(const Arguments &args, SBSearch<DB> &sbs)
    {
        // ephemeris date range is from command line or observations date range
        auto range = sbsdb::get::observations_date_range(sbs.db());
        if ((args.action != "remove") && (!args.start_date || !args.stop_date) && (!range.first || !range.second))
            throw EphemerisError("Observations database is empty: --start and --stop are required.");

        const Date start_date = args.start_date ? args.start_date.value() : range.first.value();
        const Date stop_date = args.stop_date ? args.stop_date.value() : range.second.value();

        MovingTarget target = sbsdb::get::moving_target(sbs.db(), args.target, !args.major_body);
        if (!target.moving_target_id())
        {
            sbsdb::add::moving_target(sbs.db(), target);
            cout << "Added moving target " << target.designation()
                 << " to the database with ID " << target.moving_target_id().value()
                 << "." << endl;
        }

        int count;
        if (!args.file.empty())
            count = add_from_file(args.file, target, sbs);
        else
        {
            count = add_from_horizons(target,
                                      start_date,
                                      stop_date,
                                      args.time_step,
                                      args.cache,
                                      sbs);
        }

        if (count == 0)
            throw std::runtime_error("Empty ephemeris data.");
        else
            cout << "Added " << count << " ephemeris epochs.\n\n";
    }

    template <typename DB>
    int add_from_file(const string &file, MovingTarget &target, SBSearch<DB> &sbs)
    {
        cout << "Reading ephemeris from file " << file << ".\n";
        string table = read_file(file);
        Ephemeris eph = Ephemeris(target, Horizons::parse(table));
        sbs.add_ephemeris(eph);
        return eph.num_vertices();
    }

    template <typename DB>
    int add_from_horizons(const MovingTarget &target,
                          const Date &start_date,
                          const Date &stop_date,
                          const string &time_step,
                          bool cache,
                          SBSearch<DB> &sbs)
    {
        cout << "Fetching ephemeris from Horizons for " << target.designation() << " from "
             << start_date.iso() << " to " << stop_date.iso()
             << " with step size " << time_step
             << "." << endl;

        Ephemeris eph(target, get_from_horizons(target, "500@399", start_date, stop_date, time_step, cache));

        sbs.add_ephemeris(eph);
        return eph.num_vertices();
    }

    template void add(const Arguments &, SBSearch<sbsdb::Postgresql> &);
    template int add_from_file(const string &, MovingTarget &, SBSearch<sbsdb::Postgresql> &);
    template int add_from_horizons(const MovingTarget &, const Date &, const Date &,
                                   const string &, bool, SBSearch<sbsdb::Postgresql> &);
}