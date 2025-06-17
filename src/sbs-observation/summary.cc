#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "arguments.h"
#include "summary.h"
#include "../date.h"
#include "../logging.h"
#include "../sbsearch.h"
#include "../sbsdb/count.h"
#include "../sbsdb/get.h"
#include "../sbsdb/postgresql.h"

using std::cerr;
using std::cout;
using std::string;
using std::vector;

namespace sbsearch::sbs_observation
{
    // generate a summary of observation coverage over the date range
    template <typename DB>
    void summary(const Arguments &args, SBSearch<DB> &sbs)
    {
        vector<string> sources(args.sources);
        if (sources.empty())
            sources = sbsdb::get::sources(sbs.db());

        auto range = sbsdb::get::observations_date_range(sbs.db());
        if (!range.first)
        {
            cout << "No observations to summarize.\n";
            exit(0);
        }

        const double mjd_start = args.start_date ? args.start_date.value().mjd() : range.first.value();
        const double mjd_stop = args.stop_date ? args.stop_date.value().mjd() : range.second.value();

        if (mjd_start >= mjd_stop)
        {
            cout << "Start date is after stop date.  No observations to summarize.\n";
            exit(0);
        }

        // set up histogram parameters
        const int n_bins = 100;
        const double step = (mjd_stop - mjd_start) / n_bins;

        cout << "Summarizing observation coverage over the date range "
             << Date(mjd_start).iso() << " to " << Date(mjd_stop).iso()
             << ", " << step << " day step size.\n\n";

        ProgressPercent progress(sources.size() * n_bins);
        for (const string &source : sources)
        {
            vector<int> count(n_bins, 0);
            for (int bin = 0; bin < (n_bins - 1); bin++)
            {
                count[bin] = sbsdb::count::observations(sbs.db(),
                                                        source,
                                                        mjd_start + bin * step,
                                                        mjd_start + (bin + 1) * step);
                progress.update();
                cerr << "    \r";
                progress.status(false);
            }
            cerr << "\r            \n";

            string h(n_bins, '-');
            std::transform(count.begin(), count.end(), h.begin(), [](auto i)
                           { return (i > 0) ? '+' : '-'; });

            cout << ((source == "") ? "all sources" : source) << "  " << std::accumulate(count.begin(), count.end(), 0) << "  " << h << "\n\n";
        }
    }

    template void summary(const Arguments &, SBSearch<sbsdb::Postgresql> &);
}