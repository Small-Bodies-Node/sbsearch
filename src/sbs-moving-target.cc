#include <iostream>
#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "config.h"
#include "date.h"
#include "logging.h"
#include "moving_target.h"
#include "cli.h"
#include "sbsdb.h"
#include "sbsearch.h"
#include "table.h"
#include "util/string.h"

using namespace sbsearch;
using namespace sbsearch::cli;
using namespace sbsearch::table;
using sbsearch::SBSearch;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

struct Arguments : CommonArguments
{
    string action;
    string target;
    vector<string> alternate_names;
    bool force_remove;
    std::optional<Date> start_date, stop_date;
    bool major_body;
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
        "alternate,a", value<vector<string>>(&args.alternate_names), "alternate names for the target")(
        "major-body", bool_switch(&args.major_body), "target is a major body");

    options_description remove_options("Options for remove action");
    remove_options.add_options()(
        "force", bool_switch(&args.force_remove), "do not prompt for confirmation");

    options_description summary_options("Options for summary action");
    summary_options.add_options()(
        "start", value<std::optional<Date>>(&args.start_date), "start date for summary [YYYY-MM-DD]")(
        "stop,end", value<std::optional<Date>>(&args.stop_date), "stop date for summary [YYYY-MM-DD]");

    options_description general = get_common_options((CommonArguments *)&args);

    options_description visible("");
    visible.add(add_options).add(remove_options).add(summary_options).add(general);

    options_description all("");
    all.add(visible).add(hidden);

    options_description add_action("");
    add_action.add(add_options).add(general);

    options_description list_action("");
    list_action.add(general);

    options_description remove_action("");
    remove_action.add(remove_options).add(general);

    options_description summary_action("");
    summary_action.add(summary_options).add(general);

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
                 << "<action> is one of {add, list, remove, summary}\n"
                 << visible << "\n";
        }

        if (!vm.count("action"))
            cout << "\naction is a required argument\n";

        exit(0);
    }

    validate_common_options(vm);
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

template <typename DB>
void add(const Arguments args, SBSearch<DB> &sbs)
{
    MovingTarget target{args.target, !args.major_body};
    target.add_names(args.alternate_names.begin(), args.alternate_names.end());
    sbsdb::add::moving_target(sbs.db(), target);
    cout << "Added " << target << "\n";
}

template <typename DB>
void list(const Arguments args, SBSearch<DB> &sbs)
{
    using sbsearch::util::join;

    vector<string> designations, alternates;
    for (const MovingTarget &target : sbsdb::get::all_moving_targets(sbs.db()))
    {
        designations.emplace_back(target.designation());
        auto alt = target.alternate_names();
        alternates.push_back(join<string>({alt.begin(), alt.end()}, ", "));
    }

    table::Table tab;
    tab.add(Column("designation", "%s", designations));
    tab.add(Column("alternates", "%s", alternates));
    std::cout << tab;
}

template <typename DB>
void remove(const Arguments args, SBSearch<DB> &sbs)
{
    MovingTarget target = sbsdb::get::moving_target(sbs.db(), args.target, !args.major_body);
    if (!target.moving_target_id())
        cout << args.target << " not in the database.\n";
    else
    {
        if (args.force_remove || confirm("Remove target " + target.to_string() + ", ephemeris, and found observations?"))
        {
            int count = sbsdb::remove::found(sbs.db(), target, 0, 100'000);
            cout << "Removed " << count << " found entries." << endl;

            count = sbsdb::remove::ephemeris(sbs.db(), target);
            cout << "Removed " << count << " ephemeris points." << endl;

            sbsdb::remove::moving_target(sbs.db(), target);
            cout << "Removed target " << target << "." << endl;
        }
    }
}

template <typename DB>
void summary(const Arguments args, SBSearch<DB> &sbs)
{
    // generate a summary of the ephemeris coverage of the date range
    auto range = sbsdb::get::observations_date_range(sbs.db());

    double mjd_start = args.start_date ? args.start_date.value().mjd() : range.first.value();
    double mjd_stop = args.stop_date ? args.stop_date.value().mjd() : range.second.value();

    if (mjd_start >= mjd_stop)
        mjd_stop = mjd_start + 1; // avoid rounding funniness

    vector<double> bin_edges(size_t(101));
    double step = (mjd_stop - mjd_start) / 100.;
    double mjd = mjd_start - step;
    std::generate(bin_edges.begin(), bin_edges.end(),
                  [&mjd, step]()
                  { mjd += step; return mjd; });
    bin_edges[0] = mjd_start;  // avoid rounding funniness
    bin_edges[100] = mjd_stop; // avoid rounding funniness

    cout << "Summarizing ephemeris coverage over the date range "
         << Date(mjd_start).iso() << " to " << Date(mjd_stop).iso()
         << ", " << step << " day step size.\n\n";

    auto histogram = [bin_edges](const vector<optional<double>> mjds)
    {
        vector<int> count(size_t(100), 0);
        for (auto const &mjd : mjds)
        {
            if (!mjd.has_value())
                continue;

            int i = std::upper_bound(bin_edges.begin(), bin_edges.end(), mjd.value()) - bin_edges.begin();
            if ((i > 0) && (i <= 101))
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
    for (const MovingTarget &target : sbsdb::get::all_moving_targets(sbs.db()))
    {
        string h = histogram(sbsdb::get::ephemeris(sbs.db(), target).mjd());
        cout << std::setw(16) << target.moving_target_id().value_or(-1) << "  "
             << std::setw(14) << target.designation() << "  "
             << std::setw(100) << h << "\n";
    }
}

template <typename DB>
void sbs_moving_target(int argc, char *argv[])
{
    Arguments args = get_arguments(argc, argv);

    SBSearch<DB> sbs(args.database, {args.log_file, args.log_level()});
    Logger::info() << "SBSearch moving target management tool." << endl;

    if (args.action == "add")
        add(args, sbs);
    if (args.action == "list")
        list(args, sbs);
    else if (args.action == "remove")
        remove(args, sbs);
    else if (args.action == "summary")
        summary(args, sbs);
}

int main(int argc, char *argv[])
{
    try
    {
        // get basic CLI stuff first to tell us which flavor of database to use
        string database = get_arguments(argc, argv).database;
        if (database.substr(0, 8) == "postgres")
            sbs_moving_target<sbsdb::Postgresql>(argc, argv);
        else
            throw DatabaseError("Cannot determine database type from string " + database);
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
