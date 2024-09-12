#include <iostream>
#include <string>
#include <vector>
#include <boost/json.hpp>
#include <boost/program_options.hpp>

#include "config.h"
#include "logging.h"
#include "moving_target.h"
#include "sbsearch.h"
#include "sbs-cli.h"
#include "table.h"

using namespace sbsearch;
using namespace sbsearch::cli;
using sbsearch::table::Table;
using std::cerr;
using std::cout;
using std::string;
using std::vector;
namespace json = boost::json;

struct Arguments
{
    string target;
    string input_file;
    bool small_body;

    vector<string> sources; // not yet implemented
    bool list;
    string output_filename;
    OutputFormat output_format = TableFormat;

    string database;
    string database_type;
    string log_file;
    bool verbose;
};

Arguments get_arguments(int argc, char *argv[])
{
    using namespace boost::program_options;

    Arguments args;

    options_description options("Found catalog options");
    options.add_options()(
        "target", value<string>(&args.target), "limit to this target name")(
        "input,i", value<string>(&args.input_file), "read target names from an input file")(
        "major-body", bool_switch(&args.small_body)->default_value(true), "moving target is a major body (applies to all targets in the input file)")(
        "source,s", value<vector<string>>(&args.sources), "only show results for this source data set, may be specified multiple times")(
        "list", bool_switch(&args.list), "list found observations")(
        "output,o", value<string>(&args.output_filename), "save the results to this file")(
        "format,f", value<OutputFormat>(&args.output_format), "output file format: table (default) or json");

    options_description general("General options");
    general.add_options()(
        "database,D", value<string>(&args.database)->default_value("sbsearch.db"), "SBSearch database name or file")(
        "db-type,T", value<string>(&args.database_type)->default_value("sqlite3"), "database type")(
        "log-file,L", value<string>(&args.log_file)->default_value("sbsearch.log"), "log file name")(
        "help,h", "display this help and exit")(
        "version", "output version information and exit")(
        "verbose,v", bool_switch(&args.verbose), "show debugging messages");

    options_description visible("");
    visible.add(options).add(general);

    options_description all("");
    all.add(visible);

    variables_map vm;
    boost::program_options::store(command_line_parser(argc, argv).options(all).run(), vm);
    boost::program_options::notify(vm);

    if (vm.count("version"))
    {
        cout << "SBSearch version " << SBSEARCH_VERSION << "\n";
        exit(0);
    }

    if (vm.count("help"))
    {
        cout << R"(
Usage: sbs-found [options...]

Inspect the found object catalog.
)"
             << visible << "\n";

        exit(0);
    }
    return args;
}

void list_observations(std::ostream *os, const vector<MovingTarget> targets, Arguments &args, SBSearch &sbs)
{
    Founds founds;
    for (MovingTarget target : targets)
    {
        if (target.moving_target_id() == UNDEF_MOVING_TARGET_ID)
            continue;

        founds.append(sbs.db()->get_found(target));
    }

    if (args.output_format == TableFormat)
        *os << founds;
    else // JSONFormat
    {
        json::array results = founds.as_json();
        *os << results;
    }
}

void summarize_observations(std::ostream *os, const vector<MovingTarget> targets, Arguments &args, SBSearch &sbs)
{
    vector<string> names;
    vector<double> first, last, min_rh, max_rh;
    vector<int> count;

    Founds founds;
    for (MovingTarget target : targets)
    {
        if (target.moving_target_id() == UNDEF_MOVING_TARGET_ID)
            continue;

        founds = sbs.db()->get_found(target);
        names.push_back(target.designation());
        count.push_back(founds.size());
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
                                 ->ephemeris.rh()[0]);
            max_rh.push_back(std::max_element(founds.begin(), founds.end(),
                                              [](const Found &a, const Found &b)
                                              { return a.ephemeris.rh()[0] < b.ephemeris.rh()[0]; })
                                 ->ephemeris.rh()[0]);
        }
        else
        {
            first.push_back(0);
            last.push_back(0);
            min_rh.push_back(0);
            max_rh.push_back(0);
        }
    }

    if (args.output_format == TableFormat)
    {
        Table summary;
        summary.add_column("target", "%s", names);
        summary.add_column("count", "%d", count);
        summary.add_column("first", "%.6lf", first);
        summary.add_column("last", "%.6lf", last);
        summary.add_column("min rh", "%.3f", min_rh);
        summary.add_column("max rh", "%.3f", max_rh);
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

int main(int argc, char *argv[])
{
    try
    {
        Arguments args = get_arguments(argc, argv);

        // Set log level
        int log_level = sbsearch::ERROR;
        if (args.verbose)
            log_level = sbsearch::DEBUG;

        SBSearch sbs(SBSearch::sqlite3, args.database, {args.log_file, log_level});
        Logger::info() << "SBSearch moving target found catalog tool." << std::endl;

        // Set up output stream: file or stdout
        std::ostream *os;
        std::ofstream outf;
        if (args.output_filename.empty())
            os = &cout;
        else
        {
            outf.open(args.output_filename);
            os = &outf;
        }

        // Make a list of moving targets of interest
        vector<MovingTarget> targets;
        if (!args.input_file.empty())
        {
            std::ifstream input(args.target);
            for (string name; std::getline(input, name);)
                if ((name.size() > 0) & (name[0] != '#'))
                    targets.push_back(sbs.db()->get_moving_target(name));
        }
        else if (!args.target.empty())
            targets.push_back(sbs.db()->get_moving_target(args.target));
        else
            targets = sbs.db()->get_all_moving_targets();

        // list or summarize the observations of those targets
        if (args.list)
            list_observations(os, targets, args, sbs);
        else
            summarize_observations(os, targets, args, sbs);

        if (outf.is_open())
            outf.close();
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
