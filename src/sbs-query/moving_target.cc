#include <optional>
#include <string>
#include <vector>

#include "config.h"
#include "date.h"
#include "files.h"
#include "horizons.h"
#include "ephemeris.h"
#include "logging.h"
#include "moving_target.h"
#include "sbsearch.h"
#include "sbsdb/get.h"
#include "sbsdb/postgresql.h"
#include "sbs-ephemeris/get.h"

using namespace sbsearch;
using std::cerr;
using std::cout;
using std::endl;
using std::optional;
using std::string;
using std::vector;

namespace sbsearch::sbs_query
{
    template <typename DB>
    const Founds query_moving_target(const Arguments &args, const string &designation, SBSearch<DB> &sbs)
    {
        // set up moving target
        MovingTarget target = sbsdb::get::moving_target(sbs.db(), designation, !args.major_body);

        if (!args.orbit_file.empty())
        {
            OrbitalElements orbit;
            std::ifstream inf(args.orbit_file);
            inf >> orbit;
            target.orbit(orbit);
        }

        // default is to search over all time
        const double mjd_start = args.start_date.value_or(Date(0)).mjd();
        const double mjd_stop = args.stop_date.value_or(Date(100'000)).mjd();

        Ephemeris eph;
        if (!args.eph_file.empty())
        {
            message("Reading ephemeris from file " + args.eph_file);
            eph = Ephemeris(target, Horizons::parse(read_file(args.eph_file)));
        }
        else if (args.horizons)
        {
            message("Fetching ephemeris for " + target.to_string() + " from Horizons.");

            // ephemeris date range is from command line or observations date range
            auto range = sbsdb::get::observations_date_range(sbs.db());
            if ((!args.start_date || !args.stop_date) && (!range.first || !range.second))
                throw EphemerisError("Observations database is empty: --start and --stop are required for Horizons query.");

            const Date start_date = args.start_date ? args.start_date.value() : range.first.value();
            const Date stop_date = args.stop_date ? args.stop_date.value() : range.second.value();

            // request a maximum of 1 year at a time
            auto dates = date_ranges(start_date, stop_date, 365.);
            for (auto const &[start, stop] : dates)
            {
                Ephemeris segment(target, sbs_ephemeris::get_from_horizons(target, "500@399", start_date, stop_date, args.time_step, args.cache));
                eph.append(segment);
            }
        }
        else
        {
            message("Fetching ephemeris for " + target.to_string() + " from database.");

            eph = sbsdb::get::ephemeris(sbs.db(), target, mjd_start, mjd_stop);
            if (eph.num_vertices() == 0)
                throw std::runtime_error("No ephemeris data for target found in database.");
        }

        // set up search options
        eph.mutable_options()->use_uncertainty = args.use_uncertainty;
        FindOptions find_options = {.mjd_start = mjd_start,
                                    .mjd_stop = mjd_stop,
                                    .parallax = args.parallax(),
                                    .save = args.save,
                                    .padding = args.padding,
                                    .arc_length = args.arc_length,
                                    .time_period = args.time_period,
                                    .approximate = args.approximate,
                                    .save_info = !args.info_file.empty()};

        // search
        Founds founds;
        if (args.sources.empty())
            founds = sbs.find_observations(eph, find_options);
        else
        {
            for (const string &source : args.sources)
            {
                find_options.source = source;
                founds.append(sbs.find_observations(eph, find_options));
            }
        }

        return founds;
    }

    template const Founds query_moving_target(const Arguments &, const string &, SBSearch<sbsdb::Postgresql> &);
}