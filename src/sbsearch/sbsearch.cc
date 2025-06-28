#include <string>
#include <s2/s2cell_id.h>
#include <s2/s2region_term_indexer.h>

#include "indexer.h"
#include "logging.h"
#include "sbsdb.h"
#include "sbsearch.h"

using sbsearch::sbsdb::Postgresql;
using std::string;

namespace sbsearch
{
    template <typename SBSDB>
    SBSearch<SBSDB>::SBSearch(const string &uri, const Options &options) : db_(uri)
    {
        // attempt to initialize logger
        Logger::get_logger(options.log_file).log_level(options.log_level);

        if (options.create)
            db_.setup_tables();

        indexer_ = Indexer(sbsdb::get::indexer_options(&db_));

        S2RegionTermIndexer::Options center_indexer_options;
        center_indexer_options.set_min_level(S2CellId::kMaxLevel);
        center_indexer_options.set_max_level(S2CellId::kMaxLevel);
        center_indexer_options.set_index_contains_points_only(true);
        center_indexer_ = S2RegionTermIndexer(center_indexer_options);
    };

    template SBSearch<Postgresql>::SBSearch(const string &, const Options &);
}