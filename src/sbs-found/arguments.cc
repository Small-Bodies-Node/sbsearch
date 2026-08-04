#include <iostream>
#include <stdexcept>
#include <string>
#include <boost/program_options.hpp>

#include "config.h"
#include "date.h"
#include "cli.h"
#include "sbs-found/arguments.h"

using namespace sbsearch;
using namespace sbsearch::cli;
using std::cerr;
using std::cout;
using std::endl;
using std::string;

namespace sbsearch::sbs_found
{
    Arguments get_arguments(int argc, char *argv[])
    {
        using namespace boost::program_options;

        Arguments args;

        positional_options_description positional;
        positional.add("action", 1).add("target", 1);

        options_description hidden("Hidden options");
        hidden.add_options()(
            "action", value<string>(&args.action), "action")(
            "target", value<string>(&args.target), "target");

        options_description target_options("Target options");
        target_options.add_options()(
            "input,i", bool_switch(&args.input_file), "read target names from an input file, or, when target is -, read from standard input")(
            "major-body", bool_switch(&args.major_body), "moving target is a major body (applies to all targets in the input file)");

        options_description source_options("Options for data sources");
        source_options.add_options()(
            "source,s", value<vector<string>>(&args.sources), "only find results for this source, may be specified multiple times")(
            "start", value<Date>(&args.start_date)->default_value(0),
            "start date for found data [YYYY-MM-DD or MJD]")(
            "stop,end", value<Date>(&args.stop_date)->default_value(100'000),
            "stop date for found data [YYYY-MM-DD or MJD]");

        options_description output_options("Options for output");
        output_options.add_options()(
            "output,o", value<string>(&args.output_filename), "save the results to this file")(
            "format,f", value<OutputFormat>(&args.output_format), "output file format: table (default) or json")(
            "date", value<DateFormat>(&args.date_format), "date format: mjd (default) or calendar");

        options_description remove_options("Options for remove action");
        remove_options.add_options()(
            "all", bool_switch(&args.remove_all), "remove all found results");

        options_description general = get_common_options((CommonArguments *)&args);

        options_description visible("");
        visible
            .add(target_options)
            .add(source_options)
            .add(output_options)
            .add(remove_options)
            .add(general);

        options_description all("");
        all.add(visible).add(hidden);

        options_description list_action("");
        list_action
            .add(target_options)
            .add(source_options)
            .add(output_options)
            .add(general);

        options_description remove_action("");
        remove_action
            .add(target_options)
            .add(source_options)
            .add(remove_options)
            .add(general);

        options_description summarize_action("");
        summarize_action
            .add(target_options)
            .add(source_options)
            .add(output_options)
            .add(general);

        variables_map vm;
        boost::program_options::store(
            command_line_parser(argc, argv).options(all).positional(positional).run(),
            vm);
        boost::program_options::notify(vm);

        if (vm.count("version"))
        {
            cout << "SBSearch " << SBSEARCH_VERSION << "\n";
            exit(0);
        }

        if (vm.count("help") || !vm.count("action"))
        {
            // help for a specific action?
            if (args.action == "add")
            {
                cout << "Usage: sbs-found list [target] [options...]\n\n"
                     << "List found results.\n\n"
                     << list_action << "\n";
            }
            else if (args.action == "remove")
            {
                cout << "Usage: sbs-found remove [target] [options...]\n\n"
                     << "Remove found results.\n\n"
                     << remove_action << "\n";
            }
            else if (args.action == "summarize")
            {
                cout << "Usage: sbs-found summarize [options...]\n\n"
                     << "Summarize found results.\n\n"
                     << summarize_action << "\n";
            }
            else
            {
                cout << "Usage: sbs-found <action> [target] [options...]\n\n"
                     << "Manage sbsearch found results.\n\n"
                     << "<action> is one of {list, remove, summarize}\n"
                     << visible << "\n";
            }

            if (!vm.count("action") && !vm.count("help"))
                cout << "\naction is a required argument\n";

            exit(0);
        }

        return args;
    }

}