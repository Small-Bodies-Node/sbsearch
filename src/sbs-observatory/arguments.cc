#include <iostream>
#include <string>
#include <boost/program_options.hpp>

#include "cli.h"
#include "config.h"
#include "sbs-observatory/arguments.h"

using namespace sbsearch;
using namespace sbsearch::cli;
using std::cerr;
using std::cout;
using std::string;

namespace sbsearch::sbs_observatory
{
    Arguments get_arguments(int argc, char *argv[])
    {
        using namespace boost::program_options;

        Arguments args;

        positional_options_description positional;
        positional.add("action", 1).add("name", 1).add("longitude", 1).add("rho-cos-phi", 1).add("rho-sin-phi", 1);

        options_description hidden("Hidden options");
        hidden.add_options()(
            "action", value<string>(&args.action), "target action");

        options_description add_options("Options for add action");
        add_options.add_options()(
            "name,n", value<string>(&args.observatory.name), "observatory name or IAU code for the database")(
            "longitude,l", value<double>(&args.observatory.longitude), "longitude, degrees east of Greenwich")(
            "rho-cos-phi,c", value<double>(&args.observatory.rho_cos_phi), "cosine parallax constant")(
            "rho-sin-phi,s", value<double>(&args.observatory.rho_sin_phi), "sine parallax constant");

        options_description list_options("Options for list action");
        list_options.add_options()(
            "output,o", value<string>(&args.output_filename), "save the results to this file");

        options_description general = get_common_options((CommonArguments *)&args);

        options_description visible("");
        visible.add(add_options).add(general);

        options_description all("");
        all.add(visible).add(hidden);

        options_description add_action("");
        add_action.add(add_options).add(general);

        options_description remove_action("");
        remove_action.add(general);

        variables_map vm;
        boost::program_options::store(command_line_parser(argc, argv).options(all).positional(positional).run(), vm);
        boost::program_options::notify(vm);

        if (vm.count("version"))
        {
            cout << "SBSearch version " << SBSEARCH_VERSION << "\n";
            exit(0);
        }

        if (vm.count("help") || !vm.count("action"))
        {

            // help for a specific action?
            if (args.action == "add")
            {
                cout << "Usage: sbs-observatory add <name> <longitude> <rho-cos-phi> <rho-sin-phi>\n"
                     << "Add observatory to the database.\n\n"
                     << "Coordinates are in the same format as the MPC list [1].\n"
                     << add_action << "\n"
                     << "[1] https://minorplanetcenter.net/iau/lists/ObsCodesF.html\n";
            }
            else if (args.action == "remove")
            {
                cout << "Usage: sbs-observatory remove [options...]\n"
                     << "Remove observations from the database.\n\n"
                     << remove_action << "\n";
            }
            else
            {
                cout << "Usage: sbs-observatory <action> [options...]\n\n"
                     << "Manage sbsearch observatories.\n\n"
                     << "<action> is one of {add, list, remove}\n"
                     << visible << "\n";
            }

            if (!vm.count("action"))
                cout << "\naction is a required argument\n";

            exit(0);
        }

        validate_common_options(vm);
        action_is(args.action, {"add", "list", "remove"});
        action_dependency(vm, "add", "name");
        action_dependency(vm, "add", "longitude");
        action_dependency(vm, "add", "rho-cos-phi");
        action_dependency(vm, "add", "rho-sin-phi");

        return args;
    }
}
