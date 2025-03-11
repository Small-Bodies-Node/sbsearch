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
        options.max_spatial_index_cells(get_int("SELECT value FROM configuration WHERE parameter='max_spatial_index_cells';").value());
        options.max_spatial_level(get_int("SELECT value FROM configuration WHERE parameter='max_spatial_level';").value());
        options.min_spatial_level(get_int("SELECT value FROM configuration WHERE parameter='min_spatial_level';").value());
        options.temporal_resolution(get_int("SELECT value FROM configuration WHERE parameter='temporal_resolution';").value());
        return options;
    }

    void SBSearchDatabase::update_moving_target(const MovingTarget &target)
    {
        Logger::info() << "Update moving target " << target << endl;
        remove_moving_target(target);
        MovingTarget copy(target);
        add_moving_target(copy);
        Logger::info() << target << " updated." << std::endl;
    };

    std::pair<optional<double>, optional<double>> SBSearchDatabase::ephemeris_date_range()
    {
        auto mjd_start = get_double("SELECT MIN(mjd) FROM ephemerides;");
        auto mjd_stop = get_double("SELECT MAX(mjd) FROM ephemerides;");

        return {mjd_start, mjd_stop};
    };

    void SBSearchDatabase::add_observations(Observation &observation)
    {
        Observations observations({observation});
        add_observations(observations);
        observation = observations[0];
    }
}