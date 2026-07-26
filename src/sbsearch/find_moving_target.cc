#include <algorithm>
#include <future>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>
#include <set>
#include <s2/s2polygon.h>

#include "cli.h"
#include "ephemeris.h"
#include "exceptions.h"
#include "indexer.h"
#include "logging.h"
#include "observation.h"
#include "query_info.h"
#include "queue.h"
#include "sbsdb.h"
#include "sbsearch.h"

using sbsearch::sbsdb::Postgresql;
using std::endl;
using std::optional;
using std::pair;
using std::to_string;

namespace sbsearch
{
    // Searches for observations of the ephemeris, placing unique matches in the
    // queue.
    template <class DB>
    void query_ephemeris_(const Ephemeris &ephemeris,
                          const FindOptions &options,
                          UniqueQueue<Observation, int64_t> &queue,
                          DB *db,
                          Indexer &indexer,
                          QueryInfo &query_info)
    {
        // split ephemeris into search segments
        vector<Ephemeris> segments = ephemeris.split(options.arc_length, options.time_period);
        cli::message("Ephemeris split into " + to_string(segments.size()) +
                     " segments (max " + to_string(options.arc_length) +
                     " deg, " + to_string(options.time_period) + " days per segment)");

        // search for each segment
        std::set<string> query_terms;
        ProgressPercent progress(segments.size());
        for (auto const &segment : segments)
        {
            // Account for padding and possibly parallax.
            double padding = options.padding;
            if (options.parallax)
            {
                // Increase search area by the size of the Earth at the distance
                // of the target = 8.7" / Delta = 0.145' / Delta, for Delta in au.
                auto delta = segment.delta();
                padding += 0.145 / std::max_element(delta.begin(), delta.end())->value();
            }

            // Get query terms for this segment
            vector<string> segment_query_terms = indexer.terms(Indexer::query, segment, padding);

            // Search the database given segment dates and query terms
            auto db_options = options.as_sbsearch_db_options();
            db_options.mjd_start = segment.data(0).mjd.value();
            db_options.mjd_stop = segment.data(-1).mjd.value();
            sbsdb::find::observations(db, segment_query_terms, db_options);
            Observations matches = sbsdb::find::results(db);

            // Queue up the results for detailed intersection testing.
            for (auto const &observation : matches)
                queue.put(observation, observation.observation_id().value());

            // Save the segment and matches to query_info?
            if (options.save_info)
            {
                query_info.approximate_matches(matches);
                query_info.ephemeris_segment(segment, padding, segment_query_terms);
            }

            // Report status to the user
            progress.update();
            progress.status(false);
        }
        progress.status(true);
    }

    Founds test_approximate_matches_(const Ephemeris &ephemeris,
                                     const FindOptions &options,
                                     UniqueQueue<Observation, int64_t> &queue,
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
                    observatory = observatories.at(observation.observatory());
                }
                catch (const std::out_of_range &)
                {
                    throw ObservatoryError(observation.observatory() + " not in database");
                }

                eph = eph.parallax_offset(observatory);
            }

            // only need approximate results?  done!
            if (options.approximate)
            {
                founds.append(Found(observation, eph));
                continue;
            }

            if (intersects(observation, eph, options.padding, options.mjd_start, options.mjd_stop))
                founds.data.emplace_back(observation, eph);
        }

        return std::move(founds);
    }

    template <class SBSDB>
    Founds SBSearch<SBSDB>::find_observations(const Ephemeris &ephemeris, const FindOptions &options)
    {
        options.validate();

        // reset query info
        query_info_ = QueryInfo();

        cli::message(
            "Searching for observations with ephemeris: " +
            to_string(ephemeris.as_polyline().GetLength().degrees()) + " deg, " +
            to_string(ephemeris.data(-1).mjd.value() - ephemeris.data(0).mjd.value()) + " days.");

        indexer_.mutable_options().max_spatial_query_cells(options.max_spatial_query_cells);

        // Get observatories for parallax calculations
        Observatories observatories = sbsdb::get::all_observatories(&db_);

        // queue for observations to test
        UniqueQueue<Observation, int64_t> queue;

        // run query and tests in separate thread, use the queue to share query
        // results for testing
        std::packaged_task query_task(query_ephemeris_<SBSDB>);
        std::future<void> query_result = query_task.get_future();
        std::thread query_thread(std::move(query_task),
                                 std::ref(ephemeris),
                                 std::ref(options),
                                 std::ref(queue),
                                 &db_,
                                 std::ref(indexer_),
                                 std::ref(query_info_));

        std::packaged_task testing_task(test_approximate_matches_);
        std::future<Founds> testing_result = testing_task.get_future();
        std::thread testing_thread(std::move(testing_task),
                                   std::ref(ephemeris),
                                   std::ref(options),
                                   std::ref(queue),
                                   std::ref(observatories));

        // wait for queries to complete
        query_thread.join();

        // signal that the queries are done
        queue.finish();

        // wait for testing to complete
        testing_thread.join();

        Founds founds = testing_result.get();

        // saving found items here, and not in the testing thread
        if (options.save_info)
            query_info_.matches(founds);

        char ratio[100];
        if (founds.size() == 0)
        {
            if (queue.enqueued() == 0)
                sprintf(ratio, "0:0");
            else
                sprintf(ratio, "0:1");
        }
        else
            sprintf(ratio, "%.1f:1", (float)queue.enqueued() / founds.size());

        cli::message("Matched " + to_string(founds.size()) + " of " +
                     to_string(queue.enqueued()) + " approximate matches (" +
                     ratio + ").");

        Logger::debug() << queue.total_puts() << " observations found, "
                        << queue.enqueued() << " unique observation IDs, "
                        << founds.size() << " observations found."
                        << endl;

        if (options.save)
        {
            sbsdb::add::found(&db_, founds);
            Logger::info() << founds.size() << " found observations saved to the database." << endl;
        }

        return founds;
    }

    template void query_ephemeris_(const Ephemeris &, const FindOptions &, UniqueQueue<Observation, int64_t> &, Postgresql *, Indexer &, QueryInfo &);
    template Founds SBSearch<Postgresql>::find_observations(const Ephemeris &, const FindOptions &);
}