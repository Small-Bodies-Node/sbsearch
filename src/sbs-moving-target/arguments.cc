#include <iostream>
#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "cli.h"
#include "config.h"
#include "date.h"
#include "sbs-moving-target/arguments.h"

using namespace sbsearch;
using namespace sbsearch::cli;
using sbsearch::SBSearch;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace sbsearch::sbs_moving_target
{
    Arguments get_arguments(int argc, char *argv[])
    {
        using namespace boost::program_options;

        Arguments args;

        positional_options_description positional;
        positional.add("action", 1).add("target", 1);

        options_description hidden("Hidden options");
        hidden.add_options()("action", value<string>(&args.action), "target action");
        hidden.add_options()("target", value<string>(&args.target), "target name");

        options_description add_options("Options for add action");
        add_options.add_options()(
            "alternate,a", value<vector<string>>(&args.alternate_names), "alternate names for the target")(
            "major-body", bool_switch(&args.major_body), "target is a major body");

        options_description remove_options("Options for remove action");
        remove_options.add_options()(
            "force", bool_switch(&args.force_remove), "do not prompt for confirmation");

        options_description summarize_options("Options for summarize action");
        summarize_options.add_options()(
            "start", value<std::optional<Date>>(&args.start_date), "start date for summarize [YYYY-MM-DD]")(
            "stop,end", value<std::optional<Date>>(&args.stop_date), "stop date for summarize [YYYY-MM-DD]");

        options_description general = get_common_options((CommonArguments *)&args);

        options_description visible("");
        visible.add(add_options).add(remove_options).add(summarize_options).add(general);

        options_description all("");
        all.add(visible).add(hidden);

        options_description add_action("");
        add_action.add(add_options).add(general);

        options_description list_action("");
        list_action.add(general);

        options_description remove_action("");
        remove_action.add(remove_options).add(general);

        options_description summarize_action("");
        summarize_action.add(summarize_options).add(general);

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
                cout << "Usage: sbs-moving-target add <target> [options...]\n"
                     << "Add a new moving target to the database.\n\n"
                     << "<target> is the target's primary designation\n"
                     << add_action << "\n";
            }
            else if (args.action == "list")
            {
                cout << "Usage: sbs-moving-target list [target] [options...]\n"
                     << "List all moving targets in the database.\n\n"
                     << list_action << "\n";
            }
            else if (args.action == "remove")
            {
                cout << "Usage: sbs-moving-target remove <target> [options...]\n"
                     << "Remove a moving target from the database.\n\n"
                     << "<target> is the target's name / designation\n"
                     << remove_action << "\n";
            }
            else if (args.action == "summarize")
            {
                cout << "Usage: sbs-moving-target summarize [target] [options...]\n"
                     << "Summarize the moving target database.  If no target is specified, then\n"
                     << "all targets are summarized.\n\n"
                     << summarize_action << "\n";
            }
            else
            {
                cout << "Usage: sbs-moving-target <action> [options...]\n\n"
                     << "Manage sbsearch moving targets.\n\n"
                     << "<action> is one of {add, list, remove, summarize}\n"
                     << visible << "\n";
            }

            if (!vm.count("action"))
                cout << "\naction is a required argument\n";

            exit(0);
        }

        validate_common_options(vm);
        action_is(args.action, {"add", "list", "remove", "summarize"});
        action_dependency(vm, "add", "target");
        action_dependency(vm, "remove", "target");
        action_conflicting_option(vm, "add", "start");
        action_conflicting_option(vm, "add", "stop");
        action_conflicting_option(vm, "add", "force");
        action_conflicting_option(vm, "remove", "start");
        action_conflicting_option(vm, "remove", "stop");
        action_conflicting_option(vm, "remove", "alternate");

        return args;
    }
}
