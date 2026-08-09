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
#include "progress_widgets.h"
#include "sbsearch.h"
#include "sbs-ephemeris/get.h"
#include "sbs-query/moving_target.h"
#include "sbsdb/get.h"
#include "sbsdb/postgresql.h"

using namespace sbsearch;
using std::cerr;
using std::cout;
using std::endl;
using std::optional;
using std::string;
using std::vector;

namespace sbsearch::sbs_query::moving_target
{
    template <typename DB>
    Founds query(const vector<string> &target_names,
                 const Arguments &args,
                 SBSearch<DB> &sbs,
                 std::ostream *console)
    {
        if (target_names.size() > 1)
            message::info(std::to_string(target_names.size()) + " moving targets to search.\n");

        ProgressTriangle progress;

        // default is to search over all time
        const Date search_start_date = args.start_date ? args.start_date.value() : Date(0);
        const Date search_stop_date = args.stop_date ? args.stop_date.value() : Date(100'000);

        // for ephemeris generation, limit it to what is in the database, and buffer it
        // for better interpolation with short arcs
        auto mjd_range = sbsdb::get::observations_date_range(sbs.db());
        const Date eph_start_date = std::max(search_start_date,
                                             Date(mjd_range.first.value_or(0))) -
                                    2 * EPHEMERIS_TIME_STEP;
        const Date eph_stop_date = std::min(search_stop_date,
                                            Date(mjd_range.second.value_or(100'000))) +
                                   2 * EPHEMERIS_TIME_STEP;

        FindOptions find_options = {.mjd_start = search_start_date.mjd(),
                                    .mjd_stop = search_stop_date.mjd(),
                                    .parallax = args.parallax(),
                                    .save = args.save,
                                    .padding = args.padding,
                                    .arc_length = args.arc_length,
                                    .time_period = args.time_period,
                                    .approximate = args.approximate,
                                    .save_info = !args.info_file.empty()};

        Founds founds;
        if (!args.eph_file.empty())
            founds = from_ephemeris_file(target_names.front(),
                                         args.eph_file,
                                         args.sources,
                                         args.use_uncertainty,
                                         find_options,
                                         sbs,
                                         console);
        else if (!args.orbit_file.empty())
            founds = from_orbit_file(target_names.front(),
                                     args.orbit_file,
                                     eph_start_date,
                                     eph_stop_date,
                                     args.time_step,
                                     args.cache,
                                     args.sources,
                                     find_options,
                                     sbs,
                                     console);
        else
        {
            // use database or horizons
            for (const string &name : target_names)
            {
                Founds new_founds;
                if (args.horizons)
                    new_founds = from_horizons({name},
                                               eph_start_date,
                                               eph_stop_date,
                                               args.time_step,
                                               args.cache,
                                               args.sources,
                                               args.use_uncertainty,
                                               find_options,
                                               sbs,
                                               console);

                else
                    new_founds = from_database(name,
                                               eph_start_date,
                                               eph_stop_date,
                                               args.sources,
                                               args.use_uncertainty,
                                               find_options,
                                               sbs,
                                               console);

                founds.append(new_founds);

                if (target_names.size() > 1)
                    progress.update();
            }

            if (target_names.size() == 1)
                *console << founds.size() << " observations found." << endl;
            else
            {
                progress.status();
                progress.done();
            }
        }

        return founds;
    }

    template <typename DB>
    Founds from_ephemeris(Ephemeris &eph,
                          const vector<string> &sources,
                          FindOptions &find_options,
                          SBSearch<DB> &sbs,
                          std::ostream *console)
    {
        // search
        Founds founds;
        if (sources.empty())
            founds = sbs.find_observations(eph, find_options);
        else
        {
            for (const string &source : sources)
            {
                find_options.source = source;
                founds.append(sbs.find_observations(eph, find_options));
            }
        }

        return founds;
    }

    template <typename DB>
    Founds from_database(const string &name,
                         const Date &eph_start_date,
                         const Date &eph_stop_date,
                         const vector<string> &sources,
                         const bool use_uncertainty,
                         FindOptions &find_options,
                         SBSearch<DB> &sbs,
                         std::ostream *console)
    {
        message::write("Fetching ephemeris for " + name + " from database.", *console, Logger::info());

        MovingTarget target = sbsdb::get::moving_target(sbs.db(), name);
        Ephemeris eph = sbsdb::get::ephemeris(sbs.db(), target, eph_start_date.mjd(), eph_stop_date.mjd());
        eph.mutable_options()->use_uncertainty = use_uncertainty;

        if (eph.num_vertices() == 0)
        {
            message::write("No ephemeris data for target found in database over requested time period.",
                           *console,
                           Logger::info());
            return {};
        }

        return from_ephemeris(eph, sources, find_options, sbs, console);
    }

    template <typename DB>
    Founds from_ephemeris_file(const string &name,
                               const string &eph_file,
                               const vector<string> &sources,
                               const bool use_uncertainty,
                               FindOptions &find_options,
                               SBSearch<DB> &sbs,
                               std::ostream *console)
    {
        message::write("Reading ephemeris from file " + eph_file,
                       *console, Logger::info());
        MovingTarget target(name);
        Ephemeris eph = Ephemeris(target, Horizons::parse(read_file(eph_file)));
        eph.mutable_options()->use_uncertainty = use_uncertainty;
        return from_ephemeris(eph, sources, find_options, sbs, console);
    }

    template <typename DB>
    Founds from_orbit_file(const string &name,
                           const string &orbit_file,
                           const Date &eph_start_date,
                           const Date &eph_stop_date,
                           const string &time_step,
                           bool cache,
                           const vector<string> &sources,
                           FindOptions &find_options,
                           SBSearch<DB> &sbs,
                           std::ostream *console)
    {
        message::write("Reading orbit from file " + orbit_file,
                       *console, Logger::info());

        OrbitalElements orbit;
        std::ifstream inf(orbit_file);
        inf >> orbit;
        MovingTarget target(name, orbit);

        return from_horizons(target, eph_start_date, eph_stop_date, time_step, cache, sources, false, find_options, sbs, console);
    }

    template <typename DB>
    Founds from_horizons(const MovingTarget &target,
                         const Date &eph_start_date,
                         const Date &eph_stop_date,
                         const string &time_step,
                         const bool cache,
                         const vector<string> &sources,
                         const bool use_uncertainty,
                         FindOptions &find_options,
                         SBSearch<DB> &sbs,
                         std::ostream *console)
    {
        message::write("Fetching ephemeris for " + target.to_string() + " from Horizons.",
                       *console, Logger::info());

        Ephemeris::Data data = sbs_ephemeris::get_from_horizons(target, "500@399", eph_start_date, eph_stop_date, time_step, cache);
        Ephemeris eph(target, data);
        eph.mutable_options()->use_uncertainty = use_uncertainty;

        return from_ephemeris(eph, sources, find_options, sbs, console);
    }

    template Founds query(const vector<string> &,
                          const Arguments &,
                          SBSearch<sbsdb::Postgresql> &,
                          std::ostream *);
    template Founds from_ephemeris(Ephemeris &,
                                   const vector<string> &,
                                   FindOptions &,
                                   SBSearch<sbsdb::Postgresql> &,
                                   std::ostream *);
    template Founds from_database(const string &,
                                  const Date &,
                                  const Date &,
                                  const vector<string> &,
                                  const bool,
                                  FindOptions &,
                                  SBSearch<sbsdb::Postgresql> &,
                                  std::ostream *);
    template Founds from_ephemeris_file(const string &,
                                        const string &,
                                        const vector<string> &,
                                        const bool,
                                        FindOptions &,
                                        SBSearch<sbsdb::Postgresql> &,
                                        std::ostream *);
    template Founds from_orbit_file(const string &,
                                    const string &,
                                    const Date &,
                                    const Date &,
                                    const string &,
                                    const bool,
                                    const vector<string> &,
                                    FindOptions &,
                                    SBSearch<sbsdb::Postgresql> &,
                                    std::ostream *);
    template Founds from_horizons(const MovingTarget &,
                                  const Date &,
                                  const Date &,
                                  const string &,
                                  const bool,
                                  const vector<string> &,
                                  const bool,
                                  FindOptions &,
                                  SBSearch<sbsdb::Postgresql> &,
                                  std::ostream *);
}
