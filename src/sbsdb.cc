#include "config.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <s2/s2point.h>
#include <s2/s2cap.h>
#include <s2/s2metrics.h>
#include <s2/s2polyline.h>
#include <s2/s2polygon.h>
#include <s2/s2region_term_indexer.h>

#include "ephemeris.h"
#include "indexer.h"
#include "logging.h"
#include "observation.h"
#include "observatory.h"
#include "sbsdb.h"
#include "util.h"

using std::cout;
using std::endl;
using std::string;
using std::unique_ptr;
using std::vector;

namespace sbsearch
{
    Indexer::Options SBSearchDatabase::indexer_options()
    {
        Indexer::Options options;
        options.max_spatial_index_cells(*get_int("SELECT value FROM configuration WHERE parameter=\"max_spatial_index_cells\";"));
        options.max_spatial_level(*get_int("SELECT value FROM configuration WHERE parameter=\"max_spatial_level\";"));
        options.min_spatial_level(*get_int("SELECT value FROM configuration WHERE parameter=\"min_spatial_level\";"));
        options.temporal_resolution(*get_int("SELECT value FROM configuration WHERE parameter=\"temporal_resolution\";"));
        return options;
    }

    void SBSearchDatabase::add_observations(Observations &observations) const
    {
        execute_sql("BEGIN TRANSACTION;");
        for (Observation &observation : observations)
        {
            try
            {
                add_observation(observation);
            }
            catch (std::exception &e)
            {
                Logger::error() << "Error processing observation: " << observation << std::endl;
                throw;
            }
        }
        execute_sql("END TRANSACTION;");
    }

    void SBSearchDatabase::add_founds(const Founds &founds) const
    {
        execute_sql("BEGIN TRANSACTION;");
        for (const Found &found : founds.data)
            add_found(found);
        execute_sql("END TRANSACTION;");
    }

    void SBSearchDatabase::remove_founds(const Founds &founds) const
    {
        execute_sql("BEGIN TRANSACTION;");
        for (const Found &found : founds.data)
            remove_found(found);

        execute_sql("END TRANSACTION;");
    }
}