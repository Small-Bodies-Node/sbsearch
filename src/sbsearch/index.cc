#include <iostream>
#include <string>
#include <vector>
#include <s2/s2polygon.h>

#include "sbsearch.h"
#include "../cli.h"
#include "../indexer.h"
#include "../logging.h"
#include "../observation.h"
#include "../sbsdb/sbsdb.h"
#include "../util/polygon.h"
#include "../util/string.h"

using sbsearch::sbsdb::Postgresql;
using std::endl;
using std::string;
using std::vector;

namespace sbsearch
{
    template <typename SBSDB>
    void SBSearch<SBSDB>::reindex_database_terms(const Indexer::Options &options)
    {
        sbsdb::update::indexer_options(&db_, options);
        Logger::warning() << "Database configuration has been updated." << endl;
        indexer_ = Indexer(options);

        auto n = sbsdb::count::observations(&db_, 0, 100000);
        if (n == 0)
        {
            cli::message("No observations to re-index.");
            return;
        }

        Logger::info() << "Re-indexing " << n << " observations." << endl;

        db_.drop_observations_indices();

        vector<int64_t> observation_ids;
        vector<vector<string>> observation_terms;

        const int chunk = 10000;
        observation_ids.reserve(chunk);
        observation_terms.reserve(chunk);

        ProgressPercent widget(n);
        widget.status(false);

        int64_t offset = 0;
        do
        {
            observation_ids.clear();
            observation_terms.clear();

            S2Polygon polygon;
            // get the next chunk of ids and terms
            for (auto [observation_id, fov] : sbsdb::get::all_observations_fov(&db_, chunk, offset))
            {
                observation_ids.emplace_back(observation_id);
                util::make_polygon_simple(util::make_vertices(fov), polygon);
                auto terms = indexer_.terms(Indexer::index, polygon);
                observation_terms.emplace_back(terms);
            }

            // update database terms
            sbsdb::update::observations(&db_, observation_ids, observation_terms);
            offset += observation_ids.size();
            widget += observation_ids.size();
            widget.status(false);
        } while (observation_ids.size() > 0);

        std::cout << "\nCreating observation indices." << endl;
        db_.create_observations_indices();
        Logger::info() << "\n\nRe-indexed " << widget.count() << " observations." << endl;
    }

    template <typename SBSDB>
    void SBSearch<SBSDB>::index_observations(Observations &observations)
    {
        Logger::debug() << "Indexing " << observations.size() << " observations." << endl;
        for (auto observation = observations.begin(); observation < observations.end(); observation++)
        {
            if (observation->terms().size() == 0)
                observation->terms(indexer_.terms(Indexer::index, *observation));

            if (!observation->center())
            {
                vector<S2Point> points = util::make_vertices(observation->fov());
                const S2Point point = (points[0] / 4 + points[1] / 4 + points[2] / 4 + points[3] / 4).Normalize();
                observation->center(center_indexer_.GetIndexTerms(point, "")[0]);
            }
        }
    }

    template void SBSearch<Postgresql>::reindex_database_terms(const Indexer::Options &);
    template void SBSearch<Postgresql>::index_observations(Observations &);
}