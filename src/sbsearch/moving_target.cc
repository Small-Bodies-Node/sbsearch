#include <algorithm>
#include <future>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>
#include <set>

#include "cli.h"
#include "exceptions.h"
#include "indexer.h"
#include "logging.h"
#include "observation.h"
#include "polyline.h"
#include "progress_widgets.h"
#include "query_info.h"
#include "queue.h"
#include "sbsdb.h"
#include "sbsearch.h"
#include "ephemeris/ephemeris.h"
#include "ephemeris/parallax_offset.h"
#include "ephemeris/split.h"
#include "ephemeris/safe_sampler.h"
#include "ephemeris/subsample.h"

using sbsearch::ephemeris::Ephemeris;
using sbsearch::ephemeris::SafeSampler;
using sbsearch::sbsdb::Postgresql;
using std::endl;
using std::optional;
using std::pair;
using std::to_string;

namespace sbsearch
{
    std::string make_summary_(int found, int64_t matched, int64_t unique);

    Founds test_approximate_matches_(const SafeSampler &ephemeris,
                                     const SearchOptions &options,
                                     UniqueQueue<int64_t, Observation> &queue,
                                     const Observatories &observatories);

    template <class SBSDB>
    Founds SBSearch<SBSDB>::search_observations(const Ephemeris &ephemeris,
                                                const SearchOptions &options,
                                                std::ostream &console)
    {
        // reset query info
        query_info_ = QueryInfo();

        cli::message::debug(
            "Searching for observations with ephemeris: " +
            to_string(make_polyline(ephemeris).GetLength().degrees()) + " deg, " +
            to_string(ephemeris.data().back().mjd - ephemeris.data().front().mjd) + " days.");

        // split ephemeris into search segments and add to the queue
        Queue<Ephemeris> ephemeris_queue;
        vector<Ephemeris> segments = ephemeris::split(ephemeris,
                                                      options.arc_length,
                                                      options.time_period);
        for (auto const &segment : segments)
            ephemeris_queue.put(segment);

        ephemeris_queue.finish();

        cli::message::debug("Ephemeris split into " + to_string(segments.size()) +
                            " segments (max " + to_string(options.arc_length) +
                            " deg, " + to_string(options.time_period) + " days per segment)");

        return search_observations(ephemeris_queue, options, console);
    }

    template <class SBSDB>
    Founds SBSearch<SBSDB>::search_observations(Queue<Ephemeris> &ephemeris_queue,
                                                const SearchOptions &options,
                                                std::ostream &console)
    {
        options.validate();
        indexer_.mutable_options().max_spatial_query_cells(options.max_spatial_query_cells);

        // Get observatories for parallax calculations
        Observatories observatories = sbsdb::get::all_observatories(&db_);

        // queue for observations to test
        UniqueQueue<int64_t, Observation> observations_queue;

        // Ephemeris for observations testing
        SafeSampler ephemeris;

        // run tests in separate threads
        vector<std::future<Founds>> results;
        vector<std::thread> workers;
        vector<std::packaged_task<Founds(const SafeSampler &,
                                         const SearchOptions &,
                                         UniqueQueue<int64_t, Observation> &,
                                         const Observatories &)>>
            tasks;

        for (unsigned char i = 0; i < threads(); i++)
        {
            tasks.emplace_back(test_approximate_matches_);
            results.push_back(tasks.back().get_future());
        }

        for (auto &task : tasks)
            workers.emplace_back(std::move(task),
                                 std::ref(ephemeris),
                                 options,
                                 std::ref(observations_queue),
                                 std::ref(observatories));

        // Report progress
        ProgressFraction progress(ephemeris_queue.enqueued(), console);
        progress.prefix("Ephemeris segments tested/queued: ");
        progress.status();

        // search for each segment
        std::set<string> query_terms;
        bool running = true;

        while (running)
        {
            std::optional<Ephemeris> next = ephemeris_queue.next();

            // No item... is it finished?  Done!  Else, repeat!
            if (!next)
            {
                if (ephemeris_queue.finished())
                    running = false;
                continue;
            }

            Ephemeris segment = next.value();

            // append it to the testing ephemeris object
            ephemeris.append(segment);

            // Skip this segment if it doesn't overlap with the requested time period.
            if (segment.data(-1).mjd < options.mjd_start ||
                segment.data(0).mjd > options.mjd_stop)
                continue;

            // Areal search accounts for padding and possibly parallax.
            double padding = options.padding;
            if (options.parallax)
            {
                // Increase search area by the size of the Earth at the distance
                // of the target = 8.7" / Delta = 0.145' / Delta, for Delta in au.
                auto delta = segment.delta();
                padding += 0.145 / *std::min_element(delta.begin(), delta.end());
            }

            // Get query terms for this segment
            vector<string> segment_query_terms = indexer_.terms(Indexer::query,
                                                                segment,
                                                                options.use_ephemeris_uncertainty,
                                                                padding);

            // Search the database given segment dates and query terms
            auto db_options = options.as_sbsearch_db_options();
            db_options.mjd_start = std::max(segment.data(0).mjd, options.mjd_start);
            db_options.mjd_stop = std::min(segment.data(-1).mjd, options.mjd_stop);

            sbsdb::search::observations(db(), segment_query_terms, db_options);
            Observations matches = sbsdb::search::results(db());

            // Queue up the results for detailed intersection testing.
            for (auto const &observation : matches)
                observations_queue.put(observation.observation_id().value(),
                                       observation);

            // Save the segment and matches to query_info?
            if (options.save_info)
            {
                query_info_.approximate_matches(matches);
                query_info_.ephemeris_segment(segment,
                                              options.use_ephemeris_uncertainty,
                                              padding,
                                              segment_query_terms);
            }

            // only report every N queued+tested segments
            progress.total_count(ephemeris_queue.enqueued());
            progress.update();
            if ((progress.total_count() + progress.count()) % 20 == 0)
                progress.status();
        }

        // Final progress report, but only if there wasn't one just reported
        if ((progress.total_count() + progress.count()) % 20 != 0)
            progress.status();

        // signal that the queries are done
        observations_queue.finish();

        // wait for testing to complete and get the results
        Founds founds;
        for (unsigned char i = 0; i < threads(); i++)
        {
            workers[i].join();
            founds.append(results[i].get());
        }

        // save found items to query_info?
        if (options.save_info)
            query_info_.matches(founds);

        cli::message::write(make_summary_(founds.size(),
                                          observations_queue.total_puts(),
                                          observations_queue.enqueued()),
                            console,
                            Logger::info());

        if (options.save)
        {
            sbsdb::add::found(&db_, founds);
            Logger::info() << founds.size() << " found observations saved to the database." << endl;
        }
        return founds;
    }

    // generate a nice summary message
    std::string make_summary_(int found, int64_t matched, int64_t unique)
    {
        std::stringstream stream;
        stream << matched << " nearby observations ("
               << unique << " unique), "
               << found << " matches (";

        if (found == 0)
        {
            if (unique == 0)
                stream << "0:0).";
            else
                stream << "0:1).";
        }
        else
            stream << std::setprecision(3) << (double)unique / found << ":1).";

        return stream.str();
    }

    Founds test_approximate_matches_(const SafeSampler &ephemeris,
                                     const SearchOptions &options,
                                     UniqueQueue<int64_t, Observation> &queue,
                                     const Observatories &observatories)
    {
        // the results
        Founds founds;

        bool running = true;
        while (running)
        {
            auto next = queue.next();

            // No item... is it finished?  Done!  Else, repeat!
            if (!next)
            {
                if (queue.finished())
                    running = false;
                continue;
            }
            Observation observation = next.value();

            // check for detailed intersection between ephemeris and candidates,
            // accounting for parallax
            //
            // First, get the ephemeris at the time of this observation.
            Ephemeris eph;
            try
            {
                eph = ephemeris.subsample(observation.mjd_start(), observation.mjd_stop());
            }
            catch (const std::runtime_error &)
            {
                continue;
            }

            // Account for parallax?  Then offset the ephemeris.
            if (options.parallax)
            {
                Observatory observatory;
                try
                {
                    observatory = observatories.at(string(observation.observatory()));
                }
                catch (const std::out_of_range &)
                {
                    throw ObservatoryError(string(observation.observatory()) + " not in database");
                }

                eph = ephemeris::parallax_offset(eph, observatory);
            }

            // only need approximate results?  done!
            if (options.approximate)
            {
                founds.data.emplace_back(observation, eph);
                continue;
            }

            if (intersects(observation,
                           eph,
                           options.use_ephemeris_uncertainty,
                           options.padding,
                           options.mjd_start,
                           options.mjd_stop))
                founds.data.emplace_back(observation, eph);
        }

        return std::move(founds);
    }

    template Founds SBSearch<Postgresql>::search_observations(const Ephemeris &,
                                                              const SearchOptions &,
                                                              std::ostream &);
    template Founds SBSearch<Postgresql>::search_observations(Queue<Ephemeris> &,
                                                              const SearchOptions &,
                                                              std::ostream &);
}