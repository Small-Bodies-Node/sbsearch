#include <iostream>
#include <optional>
#include <string>
#include <vector>
#include <s2/s2latlng.h>
#include <s2/s2point.h>

#include "config.h"
#include "date.h"
#include "observation.h"
#include "sbsearch.h"
#include "sbsdb/postgresql.h"
#include "sbs-query/arguments.h"
#include "sbs-query/fixed_target.h"

namespace sbsearch::sbs_query
{
    template <typename DB>
    const Observations query_fixed_target(const Arguments &args, const string &coordinates, SBSearch<DB> &sbs)
    {
        // convert target coordinates into S2Point
        const int delimiter = coordinates.find_first_of(", ");
        const double ra = std::stod(coordinates.substr(0, delimiter));
        const double dec = std::stod(coordinates.substr(delimiter + 1));
        S2Point point = S2LatLng::FromDegrees(dec, ra).Normalized().ToPoint();

        // default is to search over all time
        const double mjd_start = args.start_date.value_or(Date(0)).mjd();
        const double mjd_stop = args.stop_date.value_or(Date(100000)).mjd();

        // set options and search
        FindOptions find_options = {.mjd_start = mjd_start,
                                    .mjd_stop = mjd_stop,
                                    .padding = args.padding,
                                    .approximate = args.approximate};
        if (args.padding > 0)
            find_options.intersection_type = args.intersection_type;

        Observations observations;
        observations = sbs.find_observations(point, find_options);

        return observations;
    }

    template const Observations query_fixed_target(const Arguments &, const string &, SBSearch<sbsdb::Postgresql> &);
}