#include <istream>
#include <iostream>
#include <string>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "config.h"
#include "logging.h"
#include "observation.h"
#include "sbs-cli.h"
#include "sbsearch.h"
#include "util.h"

using sbsearch::Logger;
using sbsearch::Observation;
using sbsearch::Observations;
using sbsearch::SBSearch;
using std::cerr;
using std::cout;
using std::string;
using std::vector;

struct Arguments
{
    string action;
    string file;
    bool noop;

    bool force_remove;
    vector<string> sources;
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
    positional.add("action", 1).add("file", 2);

    options_description hidden("Hidden options");
    hidden.add_options()("action", value<string>(&args.action), "target action")(
        "file", value<string>(&args.file), "read data from this JSON-formatted file");

    options_description add_options("Options for add action");
    add_options.add_options()(
        "format-help", "display help on JSON file format and exit")(
        ",n", bool_switch(&args.noop), "no-op mode, parse the file, but do not add to the database");

    options_description remove_options("Options for remove action");
    remove_options.add_options()(
        "force,f", bool_switch(&args.force_remove), "do not prompt for confirmation");

    options_description source_options("Options for data sources");
    source_options.add_options()(
        "source,s", value<vector<string>>(&args.sources), "only remove or summarize this source, may be specified multiple times")(
        "start", value<Date>(&args.start_date),
        "start date for remove or summary [YYYY-MM-DD]")(
        "stop,end", value<Date>(&args.stop_date),
        "stop date for remove or summary [YYYY-MM-DD]");

    options_description general("General options");
    general.add_options()(
        "database,D", value<string>(&args.database)->default_value("sbsearch.db"), "SBSearch database name or file")(
        "db-type,T", value<string>(&args.database_type)->default_value("sqlite3"), "database type")(
        "log-file,L", value<string>(&args.log_file)->default_value("sbsearch.log"), "log file name")(
        "help", "display this help and exit")(
        "version", "output version information and exit")(
        "verbose,v", bool_switch(&args.verbose), "show debugging messages");

    options_description visible("");
    visible.add(add_options).add(remove_options).add(source_options).add(general);

    options_description all("");
    all.add(visible).add(hidden);

    options_description add_action("");
    add_action.add(add_options).add(general);

    options_description remove_action("");
    remove_action.add(remove_options).add(source_options).add(general);

    options_description summary_action("");
    summary_action.add(source_options).add(general);

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
            cout << "Usage: sbs-observation add <filename> [options...]\n"
                 << "Add observations to the database.\n\n"
                 << "<filename> contains JSON-formatted data\n"
                 << add_action << "\n";
        }
        else if (args.action == "remove")
        {
            cout << "Usage: sbs-observation remove [options...]\n"
                 << "Remove observations from the database.\n\n"
                 << remove_action << "\n";
        }
        else if (args.action == "summary")
        {
            cout << "Usage: sbs-observation summary [options...]\n"
                 << "Summarize the observation database.\n\n"
                 << summary_action << "\n";
        }
        else
        {
            cout << "Usage: sbs-observation <action> [options...]\n\n"
                 << "Manage sbsearch observations.\n\n"
                 << "<action> is one of {add, remove, summary}\n"
                 << visible << "\n";
        }

        if (!vm.count("action"))
            cout << "\naction is a required argument\n";

        exit(0);
    }

    if (vm.count("format-help"))
    {
        cout << R"(
The JSON file format is a list of objects:

  [
    {
      "source": "Big Survey Project",
      "observatory": "I41",
      "product_id": "unique product ID",
      "mjd_start": 60000.00,
      "mjd_stop": 60000.01,
      "fov": "0:0, 1:0, 1:1, 0:1",
      "observation_id": 1
    },
    {
      ... the next observation
    },
    ... and more
  ]

Notes:
* "fov" is comma-separated RA:Dec pairs in units of degrees.
* "observation_id" is optional, but if included and it matches a record
  in the database then the database entry is updated.

)";
        exit(0);
    }

    if ((args.action == "add") & !vm.count("file"))
        throw std::logic_error("add action requires --file");

    return args;
}

// add observations from a file
void add(const Arguments &args, SBSearch &sbs, std::istream &input)
{
    Observations observations;

    // based on a conversation with ChatGPT, 2023 March 23 version.
    boost::property_tree::ptree root;
    const int batch_size = 10000;
    observations.reserve(batch_size);
    sbsearch::ProgressTriangle progress;
    int ready = 0;
    while (input)
    {
        boost::property_tree::ptree chunk;
        std::stringstream sstr;
        char buffer[4096];

        // read some data from the input stream and write it to a string stream
        input.read(buffer, sizeof(buffer));
        sstr.write(buffer, input.gcount());

        boost::property_tree::read_json(sstr, chunk);
        root.insert(std::end(root), std::begin(chunk), std::end(chunk));
        ready += chunk.size();

        if ((ready >= batch_size) | input.eof())
        {
            auto it = std::begin(root);
            std::advance(it, progress.count());

            for (; ready > 0; ++it, --ready, ++progress)
            {
                Observation obs(
                    it->second.get<string>("source"),
                    it->second.get<string>("observatory"),
                    it->second.get<string>("product_id"),
                    it->second.get<double>("mjd_start"),
                    it->second.get<double>("mjd_stop"),
                    it->second.get<string>("fov"),
                    "",
                    it->second.get("observation_id", UNDEFINED_OBSID));
                observations.push_back(obs);
            }

            if (!args.noop)
                sbs.add_observations(observations);

            observations.clear();
            ready = 0;
        }
    }

    if (args.noop)
        cout << "Processed ";
    else
        cout << "Added ";

    cout << progress.count() << " observations.\n";
}

// remove observations by source and/or date range
void remove(const Arguments &args, SBSearch &sbs)
{
    double mjd_start = (args.start_date.mjd == -1) ? 0.0 : args.start_date.mjd;
    double mjd_stop = (args.stop_date.mjd == -1) ? 80000.0 : args.stop_date.mjd;
    int64 count = 0;

    if (args.sources.empty())
    {
        count = sbs.db()->count_observations(mjd_start, mjd_stop);

        bool remove = args.force_remove;
        if (!args.force_remove)
            remove = confirm("\nRemove " + std::to_string(count) + " observations?");

        if (remove)
        {
            cout << "Removing observations...";
            sbs.db()->remove_observations(mjd_start, mjd_stop);
            cout << "done\n\n";
        }
        else
            cout << "Cancelled.\n\n";
    }
    else
    {
        for (const string &source : args.sources)
            count += sbs.db()->count_observations(source, mjd_start, mjd_stop);

        if (args.force_remove | confirm("Remove " + std::to_string(count) + " observations from " + sbsearch::join(args.sources, ", ") + "?"))
        {
            for (const string &source : args.sources)
            {
                cout << "Removing observations from " << source << "...";
                sbs.db()->remove_observations(source, mjd_start, mjd_stop);
                cout << "done\n";
            }
        }
        else
            cout << "Cancelled.\n";
    }
}

// generate a summary of observation coverage over the date range
void summary(const Arguments &args, SBSearch &sbs)
{
    vector<string> sources(args.sources);
    if (sources.empty())
        sources = sbs.db()->get_sources();

    std::pair<double *, double *> range = sbs.db()->observation_date_range();
    if (range.first == nullptr)
    {
        cout << "No observations to summarize.\n";
        exit(0);
    }

    double mjd_start = (args.start_date.mjd == -1) ? *range.first : args.start_date.mjd;
    double mjd_stop = (args.stop_date.mjd == -1) ? *range.second : args.stop_date.mjd;

    if (mjd_start >= mjd_stop)
        mjd_stop = mjd_start + 1; // avoid rounding funniness

    // set up histogram parameters
    const size_t n_bins = 100;
    const double step = (mjd_stop - mjd_start) / n_bins;

    cout << "Summarizing observation coverage over the date range "
         << sbsearch::mjd2cal(mjd_start) << " to " << sbsearch::mjd2cal(mjd_stop)
         << ", " << step << " day step size.\n\n";

    for (const string &source : sources)
    {
        vector<int> count(n_bins, 0);
        int offset = 0;
        while (true)
        {
            Observations observations;
            if (source == "")
                observations = sbs.db()->find_observations(mjd_start, mjd_stop, 1000, offset);
            else
                observations = sbs.db()->find_observations(source, mjd_start, mjd_stop, 1000, offset);

            if (observations.size() == 0)
                break;

            for (const Observation &observation : observations)
                count[static_cast<int>((observation.mjd_start() - mjd_start) / step)]++;

            offset += observations.size();
        }
        string h(n_bins, '-');
        std::transform(count.begin(), count.end(), h.begin(), [](auto i)
                       { return (i > 0) ? '+' : '-'; });

        cout << ((source == "") ? "all sources" : source) << "  " << std::accumulate(count.begin(), count.end(), 0) << "  " << h << "\n\n";
    }
}

int main(int argc, char *argv[])
{
    try
    {
        Arguments args = get_arguments(argc, argv);
        cout << "SBSearch observation management tool.\n\n";

        SBSearch sbs(SBSearch::sqlite3, args.database, args.log_file);

        // Set log level
        if (args.verbose)
            Logger::get_logger().log_level(sbsearch::DEBUG);
        else
            Logger::get_logger().log_level(sbsearch::INFO);

        if (args.action == "add") // add data to database
        {
            if (args.file == "-")
            {
                cout << "Reading observations from stdin.\n";
                add(args, sbs, std::cin);
            }
            else
            {
                cout << "Reading observations from " + args.file + ".\n";
                std::ifstream input(args.file);
                if (!input)
                    throw std::runtime_error("Error opening file: " + args.file);
                add(args, sbs, input);
                input.close();
            }
        }
        else if (args.action == "remove")
            remove(args, sbs);
        else if (args.action == "summary")
            summary(args, sbs);
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
