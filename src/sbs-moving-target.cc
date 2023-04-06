#include <iostream>
#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "config.h"
#include "logging.h"
#include "moving_target.h"
#include "sbs-cli.h"
#include "sbsearch.h"
#include "util.h"

using sbsearch::Logger;
using sbsearch::MovingTarget;
using sbsearch::SBSearch;
using std::cerr;
using std::cout;
using std::string;
using std::vector;

struct Arguments
{
    string action;
    string target;
    vector<string> alternate_names;
    bool force_remove;
    Date start_date, stop_date;
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
    positional.add("action", 1).add("target", 1);

    options_description hidden("Hidden options");
    hidden.add_options()("action", value<string>(&args.action), "target action");
    hidden.add_options()("target", value<string>(&args.target), "target name");

    options_description add_options("Options for add action");
    add_options.add_options()(
        "alternate,a", value<vector<string>>(&args.alternate_names), "alternate names for the target");

    options_description remove_options("Options for remove action");
    remove_options.add_options()(
        "force,f", bool_switch(&args.force_remove), "do not prompt for confirmation");

    options_description summary_options("Options for summary action");
    summary_options.add_options()(
        "start", value<Date>(&args.start_date),
        "start date for summary [YYYY-MM-DD]")(
        "stop,end", value<Date>(&args.stop_date),
        "stop date for summary [YYYY-MM-DD]");

    options_description general("General options");
    general.add_options()(
        "database,D", value<string>(&args.database)->default_value("sbsearch.db"), "SBSearch database name or file")(
        "db-type,T", value<string>(&args.database_type)->default_value("sqlite3"), "database type")(
        "log-file,L", value<string>(&args.log_file)->default_value("sbsearch.log"), "log file name")(
        "help", "display this help and exit")(
        "version", "output version information and exit")(
        "verbose,v", bool_switch(&args.verbose), "show debugging messages");

    options_description visible("");
    visible.add(add_options).add(remove_options).add(summary_options).add(general);

    options_description all("");
    all.add(visible).add(hidden);

    options_description add_action("");
    add_action.add(add_options).add(general);

    options_description remove_action("");
    remove_action.add(remove_options).add(general);

    options_description summary_action("");
    summary_action.add(summary_options).add(general);

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
            cout << "Usage: sbs-moving-target add <target> [options...]\n"
                 << "Add a new moving target to the database.\n\n"
                 << "<target> is the target's primary designation\n"
                 << add_action << "\n";
        }
        else if (args.action == "remove")
        {
            cout << "Usage: sbs-moving-target remove <target> [options...]\n"
                 << "Remove a moving target from the database.\n\n"
                 << "<target> is the target's name / designation\n"
                 << remove_action << "\n";
        }
        else if (args.action == "summary")
        {
            cout << "Usage: sbs-moving-target summary [target] [options...]\n"
                 << "Summarize the moving target database.  If no target is specified, then\n"
                 << "all targets are summarized.\n\n"
                 << summary_action << "\n";
        }
        else
        {
            cout << "Usage: sbs-moving-target <action> [options...]\n\n"
                 << "Manage sbsearch moving targets.\n\n"
                 << "<action> is one of {add, remove, summary}\n"
                 << visible << "\n";
        }

        if (!vm.count("action"))
            cout << "\naction is a required argument\n";

        exit(0);
    }

    if (((args.action == "add") | (args.action == "remove")) & !vm.count("target"))
        throw std::logic_error(args.action + " action requires a target name");

    return args;
}

int main(int argc, char *argv[])
{
    try
    {
        Arguments args = get_arguments(argc, argv);
        cout << "SBSearch ephemeris management tool.\n\n";

        SBSearch sbs(SBSearch::sqlite3, args.database, args.log_file);

        // Set log level
        if (args.verbose)
            Logger::get_logger()
                .log_level(sbsearch::DEBUG);
        else
            Logger::get_logger().log_level(sbsearch::INFO);

        if (args.action == "add") // add data to database
        {
            MovingTarget target{args.target};
            target.add_names(args.alternate_names.begin(), args.alternate_names.end());
            sbs.db()->add_moving_target(target);
            cout << "Added " << target << "\n";
        }
        else if (args.action == "remove")
        {
            MovingTarget target = sbs.db()->get_moving_target(args.target);
            if (target.moving_target_id() == UNDEF_MOVING_TARGET_ID)
                cout << args.target << " not in the database.\n";
            else
            {
                if (args.force_remove | confirm("Remove target " + sbsearch::to_string(target) + "?"))
                {
                    cout << "Removing " << target << "\n";
                    sbs.db()->remove_moving_target(target);
                }
            }
        }
        else if (args.action == "summary")
        {
            // generate a summary of the ephemeris coverage of the date range
            auto range = sbs.db()->ephemeris_date_range();
            double mjd_start = args.start_date.mjd;
            double mjd_stop = args.stop_date.mjd;

            if ((mjd_start == -1) & (range.first != nullptr))
                mjd_start = *range.first;
            if ((mjd_stop == -1) & (range.second != nullptr))
                mjd_stop = *range.second;

            if (mjd_start >= mjd_stop)
                mjd_stop = mjd_start + 1; // avoid rounding funniness

            vector<double> bin_edges(size_t(101));
            double step = (mjd_stop - mjd_start) / 100.;
            double mjd = mjd_start - step;
            std::generate(bin_edges.begin(), bin_edges.end(),
                          [&mjd, step]()
                          { mjd += step; return mjd; });
            bin_edges[100] = mjd_stop; // avoid rounding funniness

            cout << "Summarizing ephemeris coverage over the date range "
                 << sbsearch::mjd2cal(mjd_start) << " to " << sbsearch::mjd2cal(mjd_stop)
                 << ", " << step << " day step size.\n\n";

            auto histogram = [bin_edges](const vector<double> mjds)
            {
                vector<int> count(size_t(100), 0);
                for (const double &mjd : mjds)
                {
                    int i = std::upper_bound(bin_edges.begin(), bin_edges.end(), mjd) - bin_edges.begin();
                    if ((i > 0) & (i <= 101))
                        count[i - 1]++;
                }

                string h(100, '-');
                std::transform(count.begin(), count.end(), h.begin(), [](auto i)
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
            for (const MovingTarget &target : sbs.db()->get_all_moving_targets())
            {
                string h = histogram(sbs.db()->get_ephemeris(target).mjd());
                cout << std::setw(16) << target.moving_target_id() << "  "
                     << std::setw(14) << target.designation() << "  "
                     << std::setw(100) << h << "\n";
            }
        }
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
