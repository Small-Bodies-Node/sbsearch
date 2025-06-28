#include <algorithm>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
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

namespace sbsearch
{
    template <class DB>
    uint query_ephemeris_(const Ephemeris &ephemeris,
                          const FindOptions &options,
                          Queue<Observation> &queue,
                          DB *db,
                          Indexer &indexer,
                          QueryInfo &query_info)
    {
        // split ephemeris into search segments
        vector<Ephemeris> segments = ephemeris.split(options.arc_length, options.time_period);
        cli::message("Ephemeris split into " + std::to_string(segments.size()) +
                     " segments (max " + std::to_string(options.arc_length) +
                     " deg, " + std::to_string(options.time_period) + " days per segment)");

        // search for each segment
        std::set<string> query_terms;
        ProgressPercent progress(segments.size());
        uint count = 0;
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
            {
                queue.put(observation);
                count++;
            }

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
        std::cout << "\n";

        return count;
    }

    Founds test_approximate_matches_(const Ephemeris &ephemeris,
                                     const FindOptions &options,
                                     Queue<Observation> &queue,
                                     const Observatories &observatories)
    {
        // Avoid testing any one observation more than once
        std::set<int64_t> tested;

        // the results
        Founds founds;

        bool running = true;
        while (running)
        {
            auto next = queue.next();

            // Queue is finished?  Done!
            if (!next)
            {
                running = false;
                continue;
            }
            Observation observation = next.value();

            // already checked?  skip it
            int64_t obsid = observation.observation_id().value();
            if (tested.find(obsid) == tested.end())
                continue;

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

            // Save this to the list of tested observations.
            tested.insert(obsid);
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
            std::to_string(ephemeris.as_polyline().GetLength().degrees()) + " deg, " +
            std::to_string(ephemeris.data(-1).mjd.value() - ephemeris.data(0).mjd.value()) + " days.");

        indexer_.mutable_options().max_spatial_query_cells(options.max_spatial_query_cells);

        // Get observatories ofr parallax calculations
        Observatories observatories = sbsdb::get::all_observatories(&db_);

        // queue for observations to test
        Queue<Observation> queue;

        // run query and tests in separate thread, use the queue to share query
        // results for testing
        std::thread query_thread(query_ephemeris_<SBSDB>,
                                 std::ref(ephemeris),
                                 std::ref(options),
                                 std::ref(queue),
                                 &db_,
                                 std::ref(indexer_),
                                 std::ref(query_info_));
        std::thread testing_thread(test_approximate_matches_,
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

        uint count = 0;
        Founds founds;

        // saving found items here, and not in the testing thread
        if (options.save_info)
            query_info_.matches(founds);

        // Logger::debug() << count << " approximate matches." << endl;

        cli::message("Matched " + std::to_string(founds.size()) + " of " +
                     std::to_string(count) + " approximate matches.");

        if (options.save)
        {
            sbsdb::add::found(&db_, founds);
            Logger::info() << count << " found observations saved to the database." << endl;
        }

        return founds;
    }

    template uint query_ephemeris_(const Ephemeris &, const FindOptions &, Queue<Observation> &, Postgresql *, Indexer &, QueryInfo &);
    template Founds SBSearch<Postgresql>::find_observations(const Ephemeris &, const FindOptions &);
}