#include <future>
#include <iostream>
#include <string>
#include <thread>
#include <utility>
#include <vector>
#include <s2/s2cell_id.h>
#include <s2/s2polygon.h>
#include <s2/s2region_term_indexer.h>

#include "cli.h"
#include "indexer.h"
#include "logging.h"
#include "observation.h"
#include "progress_widgets.h"
#include "queue.h"
#include "sbsdb.h"
#include "sbsearch.h"
#include "polygons.h"
#include "util/string.h"

using sbsearch::sbsdb::Postgresql;
using std::endl;
using std::string;
using std::vector;

namespace sbsearch
{
    // Pick observation_id and FOV string from a queue, save the ID and index terms.
    std::pair<vector<int64_t>, vector<vector<string>>> index_fov_(
        Queue<std::pair<int64_t, string>> &queue,
        const Indexer::Options &options)
    {
        vector<int64_t> observations_ids;
        vector<vector<string>> observations_terms;
        Indexer indexer(options);

        S2Polygon polygon;
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

            // generate index terms
            auto const &[obsid, fov] = next.value();
            make_polygon_simple(util::make_vertices(fov), polygon);
            auto terms = indexer.terms(Indexer::index, polygon);
            observations_ids.push_back(obsid);
            observations_terms.push_back(std::move(terms));
        }

        return {observations_ids, observations_terms};
    }

    template <typename SBSDB>
    void SBSearch<SBSDB>::reindex_database_terms(const Indexer::Options &options)
    {
        sbsdb::update::indexer_options(&db_, options);
        Logger::warning() << "Database configuration has been updated." << endl;
        indexer_ = Indexer(options);

        auto n = sbsdb::count::observations(&db_, 0, 100000);
        if (n == 0)
        {
            cli::message::info("No observations to re-index.");
            return;
        }

        cli::message::info("Re-indexing " + std::to_string(n) + " observations.");

        db_.drop_observations_indices();

        const int chunk = 10000;

        ProgressPercent widget(n);
        widget.status(false);

        // number of threads to use for indexing
        unsigned int nthreads = std::thread::hardware_concurrency();

        int64_t offset = 0;
        unsigned int count;
        do
        {
            count = 0;

            // queue up the next chunk of ids and terms
            Queue<std::pair<int64_t, string>> queue;
            for (auto next : sbsdb::get::all_observations_fov(&db_, chunk, offset))
                queue.put(next);
            queue.finish();

            // process the queue in threads
            vector<std::packaged_task<
                std::pair<vector<int64_t>, vector<vector<string>>>(
                    Queue<std::pair<int64_t, string>> &, const Indexer::Options &)>>
                tasks;
            for (unsigned int i = 0; i < nthreads; i++)
                tasks.emplace_back(index_fov_);

            vector<std::future<std::pair<vector<int64_t>, vector<vector<string>>>>> results;
            for (auto &task : tasks)
                results.push_back(task.get_future());

            vector<std::thread> threads;
            for (auto &task : tasks)
                threads.emplace_back(std::move(task), std::ref(queue), options);

            // join threads, save results to the database
            for (unsigned int i = 0; i < nthreads; i++)
            {
                threads[i].join();
                auto [observations_ids, observations_terms] = results[i].get();
                sbsdb::update::observations(&db_, observations_ids, observations_terms);
                count += observations_ids.size();
            }

            offset += count;
            widget += count;
            widget.status(false);
        } while (count > 0);

        std::cout << "\nCreating observation indices." << endl;
        db_.create_observations_indices();
        Logger::info() << "\n\nRe-indexed " << widget.count() << " observations." << endl;
    }

    // Pick Observation from a queue, index the FOV and center, as needed.
    void add_missing_indices_(Queue<Observation *> &queue, const Indexer::Options &options)
    {
        Indexer indexer(options);

        S2RegionTermIndexer::Options center_indexer_options;
        center_indexer_options.set_min_level(S2CellId::kMaxLevel);
        center_indexer_options.set_max_level(S2CellId::kMaxLevel);
        center_indexer_options.set_index_contains_points_only(true);
        S2RegionTermIndexer center_indexer = S2RegionTermIndexer(center_indexer_options);

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

            // generate index terms
            Observation *observation = next.value();

            if (observation->terms().size() == 0)
                observation->terms(indexer.terms(Indexer::index, *observation));

            if (!observation->center())
            {
                vector<S2Point> points = util::make_vertices(observation->fov());
                const S2Point point = (points[0] / 4 + points[1] / 4 + points[2] / 4 + points[3] / 4).Normalize();
                observation->center(center_indexer.GetIndexTerms(point, "")[0]);
            }
        }
    }

    template <typename SBSDB>
    void SBSearch<SBSDB>::index_observations(Observations &observations)
    {
        Logger::debug() << "Indexing " << observations.size() << " observations." << endl;

        // queue up observations for indexing
        Queue<Observation *> queue;
        for (auto &observation : observations.data)
            queue.put(&observation);
        queue.finish();

        // number of threads to use for indexing
        unsigned int nthreads = std::thread::hardware_concurrency();

        vector<std::thread> threads;
        for (unsigned int i = 0; i < nthreads; i++)
            threads.emplace_back(add_missing_indices_, std::ref(queue), indexer_options());

        // join threads
        for (unsigned int i = 0; i < nthreads; i++)
            threads[i].join();
    }

    template void SBSearch<Postgresql>::reindex_database_terms(const Indexer::Options &);
    template void SBSearch<Postgresql>::index_observations(Observations &);
}