#include <algorithm>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include "config.h"
#include "date.h"
#include "moving_target.h"
#include "sbsdb.h"
#include "sbsearch.h"
#include "sbs-moving-target/summarize.h"

using namespace sbsearch;

using sbsearch::SBSearch;
using std::cerr;
using std::cout;
using std::endl;
using std::optional;
using std::string;
using std::vector;

namespace sbsearch::sbs_moving_target
{
    template <typename DB>
    void summarize(const optional<Date> &start_date,
                   const optional<Date> &stop_date,
                   SBSearch<DB> &sbs)
    {
        // generate a summary of the ephemeris coverage of the date range
        auto range = sbsdb::get::observations_date_range(sbs.db());

        double mjd_start = start_date ? start_date.value().mjd() : range.first.value();
        double mjd_stop = stop_date ? stop_date.value().mjd() : range.second.value();

        if (mjd_start >= mjd_stop)
            mjd_stop = mjd_start + 1; // avoid rounding funniness

        vector<double> bin_edges(size_t(101));
        double step = (mjd_stop - mjd_start) / 100.;
        double mjd = mjd_start - step;
        std::generate(bin_edges.begin(), bin_edges.end(),
                      [&mjd, step]()
                      { mjd += step; return mjd; });
        bin_edges[0] = mjd_start;  // avoid rounding funniness
        bin_edges[100] = mjd_stop; // avoid rounding funniness

        cout << "Summarizing ephemeris coverage over the date range "
             << Date(mjd_start).iso() << " to " << Date(mjd_stop).iso()
             << ", " << step << " day step size.\n\n";

        auto histogram = [bin_edges](const vector<double> mjds)
        {
            vector<int> count(size_t(100), 0);
            for (auto const &mjd : mjds)
            {
                int i = std::upper_bound(bin_edges.begin(), bin_edges.end(), mjd) - bin_edges.begin();
                if ((i > 0) && (i <= 101))
                    count[i - 1]++;
            }

            string h(100, '-');
            std::transform(count.begin(), count.end(), h.begin(),
                           [](auto i)
                           { return (i > 0) ? '+' : '-'; });
            return h;
        };

        cout
            << std::setw(18) << "moving_target_id  "
            << std::setw(16) << "designation  "
            << std::setw(100) << "coverage"
            << "\n"
            << std::setfill('-') << std::setw(16) << ""
            << "  "
            << std::setw(14) << ""
            << "  "
            << std::setw(100) << ""
            << "\n"
            << std::setfill(' ');
        for (const MovingTarget &target : sbsdb::get::all_moving_targets(sbs.db()))
        {
            string h = histogram(sbsdb::get::ephemeris(sbs.db(), target).mjd());
            cout << std::setw(16) << target.moving_target_id().value_or(-1) << "  "
                 << std::setw(14) << target.designation() << "  "
                 << std::setw(100) << h << "\n";
        }
    }

    template void summarize(const optional<Date> &, const optional<Date> &, SBSearch<sbsdb::Postgresql> &);
}
