#include <iostream>
#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "config.h"
#include "logging.h"
#include "moving_target.h"
#include "sbsearch.h"
#include "sbs-cli.h"

using namespace sbsearch;
using namespace sbsearch::cli;
using std::cerr;
using std::cout;
using std::string;
using std::vector;

struct Arguments
{
    string target;
    bool input_file;
    bool small_body;

    vector<string> sources; // not yet implemented
    string output_filename; // not yet implemented

    string database;
    string database_type;
    string log_file;
    bool verbose;
};

Arguments get_arguments(int argc, char *argv[])
{
    using namespace boost::program_options;

    Arguments args;

    positional_options_description positional;
    positional.add("target", 1);

    options_description hidden("Hidden options");
    hidden.add_options()("target", value<string>(&args.target), "target");

    options_description options("Found catalog options");
    options.add_options()(
        "input,i", bool_switch(&args.input_file), "read target names from an input file")(
        "major-body", bool_switch(&args.small_body)->default_value(true), "moving target is a major body (applies to all targets in the input file)")(
        "source,s", value<vector<string>>(&args.sources), "only show results for this source data set, may be specified multiple times")(
        "output,o", value<string>(&args.output_filename), "save the results to this file");

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
    all.add(visible).add(hidden);

    variables_map vm;
    boost::program_options::store(command_line_parser(argc, argv).options(all).positional(positional).run(), vm);
    boost::program_options::notify(vm);

    if (vm.count("version"))
    {
        cout << "SBSearch version " << SBSEARCH_VERSION << "\n";
        exit(0);
    }

    if (vm.count("help") | !vm.count("target"))
    {
        cout << R"(
Usage: sbs-found <target> [options...]

Inspect the found object catalog.

<target> is a moving target name or, with --input, a file
listing multiple targets.
)"
             << visible << "\n";

        if (!vm.count("target"))
            cout << "\ntarget is a required argument\n";

        exit(0);
    }
    return args;
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

        vector<string> names;
        if (args.input_file)
        {
            std::ifstream input(args.target);
            for (string line; std::getline(input, line);)
                if ((line.size() > 0) & (line[0] != '#'))
                    names.push_back(line);
        }
        else if (!args.target.empty())
            names.push_back(args.target);

        Founds founds;
        for (string name : names)
        {
            MovingTarget target = sbs.db()->get_moving_target(name);
            if (target.moving_target_id() == UNDEF_MOVING_TARGET_ID)
                continue;

            founds.append(sbs.db()->get_found(target));
        }

        cout << founds;
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
