#include <iostream>
#include <string>
#include <boost/program_options.hpp>

#include "config.h"
#include "logging.h"
#include "observatory.h"
#include "sbs-cli.h"
#include "sbsearch.h"
#include "table.h"
#include "util.h"

using namespace sbsearch;
using namespace sbsearch::cli;
using sbsearch::table::Table;
using std::cerr;
using std::cout;
using std::string;

struct Arguments
{
    string action;
    string name;
    Observatory observatory;

    string output_filename;

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
    positional.add("action", 1).add("name", 1).add("longitude", 1).add("rho-cos-phi", 1).add("rho-sin-phi", 1);

    options_description hidden("Hidden options");
    hidden.add_options()(
        "action", value<string>(&args.action), "target action");

    options_description add_options("Options for add action");
    add_options.add_options()(
        "name,n", value<string>(&args.name), "observatory name or IAU code for the database")(
        "longitude,l", value<double>(&args.observatory.longitude), "longitude, degrees east of Greenwich")(
        "rho-cos-phi,c", value<double>(&args.observatory.rho_cos_phi), "cosine parallax constant")(
        "rho-sin-phi,s", value<double>(&args.observatory.rho_sin_phi), "sine parallax constant");

    options_description list_options("Options for list action");
    list_options.add_options()(
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

    if (vm.count("help") | !vm.count("action"))
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

    action_dependency(vm, "add", "name");
    action_dependency(vm, "add", "longitude");
    action_dependency(vm, "add", "rho-cos-phi");
    action_dependency(vm, "add", "rho-sin-phi");

    return args;
}

void list(Arguments &args, SBSearch &sbs)
{
    Observatories observatories = sbs.db()->get_observatories();

    std::ostream *os;
    std::ofstream outf;

    if (args.output_filename.empty())
        os = &cout;
    else
    {
        outf.open(args.output_filename);
        os = &outf;
    }

    if (observatories.size() == 0)
        cout << "# No observatories in the database.\n";
    else
    {
        int N = observatories.size();
        vector<string> names(N);
        vector<double> lon(N), rho_cos_phi(N), rho_sin_phi(N);

        for (const auto &item : observatories)
        {
            names.push_back(item.first);
            lon.push_back(item.second.longitude);
            rho_cos_phi.push_back(item.second.rho_cos_phi);
            rho_sin_phi.push_back(item.second.rho_sin_phi);
        }

        Table table;
        table.add_column("name", "%s", names);
        table.add_column("longitude", "%.3lf", lon);
        table.add_column("rho cos(phi)", "%.6lf", rho_cos_phi);
        table.add_column("rho sin(phi)", "%.6lf", rho_sin_phi);
        *os << table;
    }
}

int main(int argc, char *argv[])
{
    try
    {
        Arguments args = get_arguments(argc, argv);

        // Set log level
        int log_level = sbsearch::INFO;
        if (args.verbose)
            log_level = sbsearch::DEBUG;

        SBSearch sbs(SBSearch::sqlite3, args.database, {args.log_file, log_level});
        Logger::info() << "SBSearch observatory management tool." << std::endl;

        if (args.action == "add") // add data to database
            sbs.db()->add_observatory(args.name, args.observatory);
        else if (args.action == "list")
            list(args, sbs);
        else if (args.action == "remove")
            sbs.db()->remove_observatory(args.name);
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
