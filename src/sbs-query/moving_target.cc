#include <future>
#include <optional>
#include <string>
#include <thread>
#include <vector>
#include <boost/json.hpp>

#include "config.h"
#include "date.h"
#include "files.h"
#include "horizons.h"
#include "logging.h"
#include "moving_target.h"
#include "progress_widgets.h"
#include "queue.h"
#include "ephemeris/ephemeris.h"
#include "ephemeris/split.h"
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

        // for ephemeris generation, limit it to what is in the database
        auto obs_mjd_range = sbsdb::get::observations_date_range(sbs.db());
        if (!obs_mjd_range.first.has_value())
        {
            Logger::warning() << "No observations in the database to search." << std::endl;
            return {};
        }

        const Date eph_start_date = std::max(search_start_date, Date(obs_mjd_range.first.value()));
        const Date eph_stop_date = std::min(search_stop_date, Date(obs_mjd_range.second.value()));

        SearchOptions search_options = {.mjd_start = search_start_date.mjd(),
                                        .mjd_stop = search_stop_date.mjd(),
                                        .parallax = args.parallax(),
                                        .save = args.save,
                                        .padding = args.padding,
                                        .use_ephemeris_uncertainty = args.use_uncertainty,
                                        .arc_length = args.arc_length,
                                        .time_period = args.time_period,
                                        .approximate = args.approximate,
                                        .save_info = !args.info_file.empty(),
                                        .debug = args.debug};

        Founds founds;
        if (!args.eph_file.empty())
            founds = from_ephemeris_file(target_names.front(),
                                         args.eph_file,
                                         args.sources,
                                         search_options,
                                         sbs,
                                         console);
        else if (!args.orbit_file.empty())
            founds = from_orbit_file(target_names.front(),
                                     args.orbit_file,
                                     eph_start_date,
                                     eph_stop_date,
                                     args.step_size,
                                     args.cache,
                                     args.sources,
                                     search_options,
                                     sbs,
                                     console);
        else
        {
            // use database or horizons
            for (string_view name : target_names)
            {
                Founds new_founds;
                if (args.horizons)
                    new_founds = from_horizons({name},
                                               eph_start_date,
                                               eph_stop_date,
                                               args.step_size,
                                               args.cache,
                                               args.sources,
                                               search_options,
                                               sbs,
                                               console);
                else
                    // Buffer database ephemerides to help get good interpolation data
                    new_founds = from_database(name,
                                               eph_start_date - 2 * EPHEMERIS_TIME_STEP,
                                               eph_stop_date + 2 * EPHEMERIS_TIME_STEP,
                                               args.sources,
                                               search_options,
                                               sbs,
                                               console);

                founds.append(new_founds);

                if (target_names.size() > 1)
                    progress.update();
            }

            if (target_names.size() != 1)
            {
                progress.done();
                std::cout << "Found " << founds.size() << " observations.\n";
            }
        }

        return founds;
    }

    template <typename DB>
    Founds from_ephemeris(Ephemeris &eph,
                          const vector<string> &sources,
                          SearchOptions &search_options,
                          SBSearch<DB> &sbs,
                          std::ostream *console)
    {
        // search
        Founds founds;
        if (sources.empty())
            founds = sbs.search_observations(eph, search_options, *console);
        else
        {
            for (string_view source : sources)
            {
                search_options.source = source;
                founds.append(sbs.search_observations(eph, search_options, *console));
            }
        }

        return founds;
    }

    template <typename DB>
    Founds from_database(string_view name,
                         const Date &eph_start_date,
                         const Date &eph_stop_date,
                         const vector<string> &sources,
                         SearchOptions &search_options,
                         SBSearch<DB> &sbs,
                         std::ostream *console)
    {
        message::write("Fetching ephemeris from database.", *console, Logger::info());

        MovingTarget target = sbsdb::get::moving_target(sbs.db(), name);
        Ephemeris eph = sbsdb::get::ephemeris(sbs.db(),
                                              target,
                                              eph_start_date.mjd(),
                                              eph_stop_date.mjd());

        if (eph.num_vertices() == 0)
        {
            message::write("No ephemeris data for target found in database over requested time period.",
                           *console,
                           Logger::info());
            return {};
        }

        return from_ephemeris(eph, sources, search_options, sbs, console);
    }

    template <typename DB>
    Founds from_ephemeris_file(string_view name,
                               string_view eph_file,
                               const vector<string> &sources,
                               SearchOptions &search_options,
                               SBSearch<DB> &sbs,
                               std::ostream *console)
    {
        message::write("Reading ephemeris from file " + string{eph_file},
                       *console,
                       Logger::info());
        MovingTarget target(name);
        Ephemeris eph = Ephemeris(target, Horizons::parse(read_file(eph_file)));
        return from_ephemeris(eph, sources, search_options, sbs, console);
    }

    template <typename DB>
    Founds from_orbit_file(string_view name,
                           string_view orbit_file,
                           const Date &eph_start_date,
                           const Date &eph_stop_date,
                           string_view step_size,
                           bool cache,
                           const vector<string> &sources,
                           SearchOptions &search_options,
                           SBSearch<DB> &sbs,
                           std::ostream *console)
    {
        namespace json = boost::json;

        message::write("Reading orbit from file " + string{orbit_file},
                       *console,
                       Logger::info());

        OrbitalElements orbit;

        std::ifstream inf(string{orbit_file});
        if ((orbit_file.size() > 5) && (orbit_file.substr(orbit_file.size() - 5) == ".json"))
        {
            json::value jv = json::parse(inf);
            orbit = OrbitalElements::from_json(jv.as_object());
        }
        else
            inf >> orbit;

        return from_horizons(MovingTarget(name, orbit),
                             eph_start_date,
                             eph_stop_date,
                             step_size,
                             cache,
                             sources,
                             search_options,
                             sbs,
                             console);
    }

    template <typename DB>
    Founds from_horizons(const MovingTarget &target,
                         const Date &eph_start_date,
                         const Date &eph_stop_date,
                         string_view step_size,
                         const bool cache,
                         const vector<string> &sources,
                         SearchOptions &search_options,
                         SBSearch<DB> &sbs,
                         std::ostream *console)
    {
        message::write("Fetching ephemeris from Horizons.", *console, Logger::info());

        // Set up a queue and worker to process ephemeris segments while we
        // query Horizons in this thread
        Queue<Ephemeris> queue;
        auto search = [&sbs, &queue, &console, &search_options]()
        {
            return sbs.search_observations(queue, search_options, *console);
        };

        std::packaged_task task(search);
        std::future<Founds> results = task.get_future();
        std::thread worker(std::move(task));

        // Split up the ephemeris request here
        vector<std::pair<Date, Date>> dates = date_ranges(eph_start_date, eph_stop_date, 365.);

        for (auto const &[start, stop] : dates)
        {
            Ephemeris::Data data = sbs_ephemeris::get_from_horizons(target,
                                                                    "500@399",
                                                                    start,
                                                                    stop,
                                                                    step_size,
                                                                    cache);

            // Split up the ephemeris itself here
            vector<Ephemeris> segments = ephemeris::split(Ephemeris(target, data),
                                                          search_options.arc_length,
                                                          search_options.time_period);
            for (auto const &segment : segments)
                queue.put(segment);
        }
        queue.finish();
        worker.join();

        return results.get();
    }

    template Founds query(const vector<string> &,
                          const Arguments &,
                          SBSearch<sbsdb::Postgresql> &,
                          std::ostream *);
    template Founds from_ephemeris(Ephemeris &,
                                   const vector<string> &,
                                   SearchOptions &,
                                   SBSearch<sbsdb::Postgresql> &,
                                   std::ostream *);
    template Founds from_database(string_view,
                                  const Date &,
                                  const Date &,
                                  const vector<string> &,
                                  SearchOptions &,
                                  SBSearch<sbsdb::Postgresql> &,
                                  std::ostream *);
    template Founds from_ephemeris_file(string_view,
                                        string_view,
                                        const vector<string> &,
                                        SearchOptions &,
                                        SBSearch<sbsdb::Postgresql> &,
                                        std::ostream *);
    template Founds from_orbit_file(string_view,
                                    string_view,
                                    const Date &,
                                    const Date &,
                                    string_view,
                                    const bool,
                                    const vector<string> &,
                                    SearchOptions &,
                                    SBSearch<sbsdb::Postgresql> &,
                                    std::ostream *);
    template Founds from_horizons(const MovingTarget &,
                                  const Date &,
                                  const Date &,
                                  string_view,
                                  const bool,
                                  const vector<string> &,
                                  SearchOptions &,
                                  SBSearch<sbsdb::Postgresql> &,
                                  std::ostream *);
}
