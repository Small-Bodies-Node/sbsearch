#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <boost/json.hpp>

#include "config.h"
#include "cli.h"
#include "found.h"
#include "moving_target.h"
#include "sbsdb.h"
#include "sbsearch.h"
#include "sbs-found/filter_sources.h"
#include "sbs-found/list.h"

using namespace sbsearch;
using namespace sbsearch::cli;
namespace json = boost::json;

using std::cout;
using std::string;
using std::vector;

namespace sbsearch::sbs_found
{
    template <typename DB>
    void list_found(SBSearch<DB> &sbs,
                    const double start_mjd,
                    const double stop_mjd,
                    const vector<string> &sources,
                    const string output_filename,
                    const OutputFormat output_format)
    {
        Founds founds = sbsdb::get::found(sbs.db(), start_mjd, stop_mjd);
        filter_sources(sources, founds);
        list_found(founds, output_filename, output_format);
    }

    template <typename DB>
    void list_found(SBSearch<DB> &sbs,
                    const vector<MovingTarget> &targets,
                    const double start_mjd,
                    const double stop_mjd,
                    const vector<string> &sources,
                    const string output_filename,
                    const OutputFormat output_format)
    {
        Founds founds;
        for (MovingTarget target : targets)
        {
            if (!target.moving_target_id())
                continue;

            founds.append(sbsdb::get::found(sbs.db(), target, start_mjd, stop_mjd));
        }

        filter_sources(sources, founds);
        list_found(founds, output_filename, output_format);
    }

    void list_found(const Founds &founds,
                    const string output_filename,
                    const OutputFormat output_format)
    {
        // Set up output stream: file or stdout
        std::ostream *os;
        std::ofstream outf;
        if (output_filename.empty())
            os = &cout;
        else
        {
            outf.open(output_filename);
            os = &outf;
        }

        if (output_format == TABLE)
            *os << founds;
        else
        {
            json::array results = founds.as_json();
            *os << results;
        }
    }

    template void list_found(SBSearch<sbsdb::Postgresql> &, const double, const double, const vector<string> &sources, const string, const OutputFormat);
    template void list_found(SBSearch<sbsdb::Postgresql> &, const vector<MovingTarget> &, const double, const double, const vector<string> &sources, const string, const OutputFormat);
}