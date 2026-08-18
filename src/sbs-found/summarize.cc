#include <iostream>
#include <string>
#include <vector>
#include <boost/json.hpp>
#include <boost/program_options.hpp>

#include "config.h"
#include "cli.h"
#include "date.h"
#include "logging.h"
#include "moving_target.h"
#include "sbsdb.h"
#include "sbsearch.h"
#include "table.h"
#include "sbs-found/filter_sources.h"
#include "sbs-found/summarize.h"

using namespace sbsearch;
using namespace sbsearch::table;
namespace json = boost::json;

using std::cerr;
using std::cout;
using std::string;
using std::vector;

namespace sbsearch::sbs_found
{
    template <typename DB>
    void summarize_found(SBSearch<DB> &sbs,
                         const vector<MovingTarget> &targets,
                         const double start_mjd,
                         const double stop_mjd,
                         const vector<string> &sources,
                         const string output_filename,
                         const cli::OutputFormat output_format)
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

        // summarize by target
        vector<string> names;
        vector<double> first, last, min_rh, max_rh;
        vector<int> count;

        Founds founds;
        for (MovingTarget target : targets)
        {
            if (!target.moving_target_id())
                continue;

            founds = sbsdb::get::found(sbs.db(), target, start_mjd, stop_mjd);
            names.emplace_back(target.designation());
            count.emplace_back(founds.size());
            if (founds.size() > 0)
            {
                first.push_back(std::min_element(founds.begin(), founds.end(),
                                                 [](const Found &a, const Found &b)
                                                 { return a.observation.mjd_mid() < b.observation.mjd_mid(); })
                                    ->observation.mjd_mid());
                last.push_back(std::max_element(founds.begin(), founds.end(),
                                                [](const Found &a, const Found &b)
                                                { return a.observation.mjd_mid() < b.observation.mjd_mid(); })
                                   ->observation.mjd_mid());
                min_rh.push_back(std::min_element(founds.begin(), founds.end(),
                                                  [](const Found &a, const Found &b)
                                                  { return a.ephemeris.rh()[0] < b.ephemeris.rh()[0]; })
                                     ->ephemeris.rh()[0]
                                     .value());
                max_rh.push_back(std::max_element(founds.begin(), founds.end(),
                                                  [](const Found &a, const Found &b)
                                                  { return a.ephemeris.rh()[0] < b.ephemeris.rh()[0]; })
                                     ->ephemeris.rh()[0]
                                     .value());
            }
            else
            {
                first.push_back(0);
                last.push_back(0);
                min_rh.push_back(0);
                max_rh.push_back(0);
            }
        }

        if (output_format == cli::OutputFormat::TABLE)
        {
            Table summary;
            summary.add(Column("target", "%s", names));
            summary.add(Column("count", "%d", count));
            summary.add(Column("first", "%.6lf", first));
            summary.add(Column("last", "%.6lf", last));
            summary.add(Column("min rh", "%.3f", min_rh));
            summary.add(Column("max rh", "%.3f", max_rh));
            *os << summary;
        }
        else
        {
            json::array output;
            for (int i = 0; i < names.size(); i++)
                output.emplace_back(
                    json::object{
                        {"target", names[i]},
                        {"count", count[i]},
                        {"first", first[i]},
                        {"last", last[i]},
                        {"min rh", min_rh[i]},
                        {"max rh", max_rh[i]}});
            *os << output;
        }
    }

    template void summarize_found(SBSearch<sbsdb::Postgresql> &,
                                  const vector<MovingTarget> &,
                                  const double,
                                  const double,
                                  const vector<string> &,
                                  const string,
                                  const cli::OutputFormat);
}