#include <iostream>
#include <optional>
#include <string>
#include <vector>
#include <s2/s2latlng.h>
#include <s2/s2point.h>

#include "config.h"
#include "date.h"
#include "observation.h"
#include "progress_widgets.h"
#include "sbsearch.h"
#include "sbsdb/postgresql.h"
#include "sbs-query/arguments.h"
#include "sbs-query/fixed_target.h"
#include "util/string.h"

namespace sbsearch::sbs_query::fixed_target
{
    template <typename DB>
    Observations query(const vector<string> &targets,
                       const Arguments &args,
                       SBSearch<DB> &sbs,
                       std::ostream *console)
    {
        // default is to search over all time
        const double mjd_start = args.start_date.value_or(Date(0)).mjd();
        const double mjd_stop = args.stop_date.value_or(Date(100'000)).mjd();

        // set options and search
        SearchOptions search_options = {.mjd_start = mjd_start,
                                        .mjd_stop = mjd_stop,
                                        .padding = args.padding,
                                        .approximate = args.approximate};
        if (args.padding > 0)
            search_options.intersection_type = args.intersection_type;

        ProgressTriangle progress;
        Observations observations;
        for (string_view target : targets)
        {
            observations.append(from_coordinates(target, search_options, sbs));
            if (targets.size() > 1)
                progress.update();
        }

        if (targets.size() > 1)
        {
            progress.status();
            progress.done();
        }

        return observations;
    }

    template <typename DB>
    Observations from_coordinates(string_view coordinates,
                                  const SearchOptions &search_options,
                                  SBSearch<DB> &sbs)
    {
        // convert target coordinates into S2Point
        const int delimiter = coordinates.find_first_of(", ");
        const double ra = util::svtod(coordinates.substr(0, delimiter));
        const double dec = util::svtod(coordinates.substr(delimiter + 1));
        S2Point point = S2LatLng::FromDegrees(dec, ra).Normalized().ToPoint();

        Observations observations;
        observations = sbs.search_observations(point, search_options);

        return observations;
    }

    template Observations query(const vector<string> &,
                                const Arguments &,
                                SBSearch<sbsdb::Postgresql> &,
                                std::ostream *);
    template Observations from_coordinates(string_view,
                                           const SearchOptions &,
                                           SBSearch<sbsdb::Postgresql> &);
}