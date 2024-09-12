#include <iostream>
#include <string>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/json.hpp>
#include <boost/system.hpp>

#include "config.h"
#include "logging.h"
#include "observation.h"
#include "sbs-cli.h"
#include "sbsearch.h"
#include "util.h"

using namespace sbsearch;
using namespace sbsearch::cli;
using std::cerr;
using std::cout;
using std::string;
using std::vector;

namespace json = boost::json;

struct Arguments
{
    string action;
    string file;
    int batch_size;
    bool drop_indices;
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
    hidden.add_options()(
        "action", value<string>(&args.action), "target action")(
        "file", value<string>(&args.file), "read data from this JSON-formatted file");

    options_description add_options("Options for add action");
    add_options.add_options()(
        "format-help", "display help on JSON file format and exit")(
        ",i", bool_switch(&args.drop_indices), "drop observations indices before adding add, re-build indices when done")(
        "batch-size,b", value<int>(&args.batch_size)->default_value(10000), "expect up to <n> observations per JSON object")(
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
        "help,h", "display this help and exit")(
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
            cout << "Usage: sbs-observation add <file> [options...]\n"
                 << "Add observations to the database.\n\n"
                 << "<file> contains JSON-formatted data\n"
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
  [ ... additional arrays as needed ]

Notes:
* "fov" is comma-separated RA:Dec pairs in units of degrees.
* "observation_id" is optional, but if included and it matches a record
  in the database then the database entry is updated.
* Multiple JSON arrays of observations may be included in the file.  This
  helps with memory conservation when >>10000 observations are being added.

)";
        exit(0);
    }

    action_dependency(vm, "add", "file");

    return args;
}

namespace sbsearch
{
    Observation tag_invoke(json::value_to_tag<Observation>, json::value const &jv)
    {
        json::object const &obj = jv.as_object();

        Observation obs(
            json::value_to<string>(obj.at("source")),
            json::value_to<string>(obj.at("observatory")),
            json::value_to<string>(obj.at("product_id")),
            json::value_to<double>(obj.at("mjd_start")),
            json::value_to<double>(obj.at("mjd_stop")),
            json::value_to<string>(obj.at("fov")));

        if (obj.contains("observation_id"))
            obs.observation_id(json::value_to<int64>(obj.at("observation_id")));

        return obs;
    }
}

// add observations from a stream
void add(const Arguments &args, SBSearch &sbs, std::istream &input)
{
    sbsearch::ProgressTriangle progress;

    Observations observations;
    observations.reserve(args.batch_size);

    json::stream_parser parser;
    boost::system::error_code error;
    parser.reset();

    if (args.drop_indices)
    {
        Logger::info() << "Dropping observations indices." << std::endl;
        sbs.drop_observations_indices();
    }

    size_t buffered = 0;
    while (input)
    {
        char buffer[4096];
        input.read(buffer, sizeof(buffer));

        size_t n = parser.write_some(buffer, input.gcount(), error);
        buffered += n;
        if (error)
            throw std::runtime_error("Error parsing JSON:" + error.message());

        int remainder = input.gcount() - n; // number of unconsumed characters

        // when JSON objects are ready, process them
        while (parser.done())
        {
            json::value data = parser.release();
            parser.reset();
            buffered = 0;

            Observations observations = json::value_to<Observations>(data);
            if (!args.noop)
                sbs.add_observations(observations);
            progress += observations.size();

            if (remainder)
            {
                int m = parser.write_some(buffer + n, remainder);
                remainder -= m; // update unconsumed characters
                buffered += m;
                n += m;
            }
        }
    }

    if (buffered > 0)
        // input is done, and parser was done parsing objects, but there is
        // still data in the buffer, this is an error:
        throw std::runtime_error("Processing complete, but " + std::to_string(buffered) + " bytes remain in the buffer.");

    cout << "Added " << progress.count() << " observations.\n\n";

    if (args.drop_indices)
    {
        Logger::info() << "Building observations indices." << std::endl;
        sbs.create_observations_indices();
    }
}

// remove observations by source and/or date range
void remove(const Arguments &args, SBSearch &sbs)
{
    double mjd_start = (args.start_date.mjd() == -1) ? 0.0 : args.start_date.mjd();
    double mjd_stop = (args.stop_date.mjd() == -1) ? 80000.0 : args.stop_date.mjd();
    int64 count = 0;

    if (args.sources.empty())
    {
        count = sbs.db()->count_observations(mjd_start, mjd_stop);
        if (!count)
        {
            cout << "No observations to remove.\n";
            return;
        }

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

    double mjd_start = (args.start_date.mjd() == -1) ? *range.first : args.start_date.mjd();
    double mjd_stop = (args.stop_date.mjd() == -1) ? *range.second : args.stop_date.mjd();

    if (mjd_start >= mjd_stop)
    {
        cout << "Start date is after stop date.  No observations to summarize.\n";
        exit(0);
    }

    // set up histogram parameters
    const int n_bins = 100;
    const double step = (mjd_stop - mjd_start) / n_bins;

    cout << "Summarizing observation coverage over the date range "
         << sbsearch::mjd2cal(mjd_start) << " to " << sbsearch::mjd2cal(mjd_stop)
         << ", " << step << " day step size.\n\n";

    ProgressPercent progress(sources.size() * n_bins);
    for (const string &source : sources)
    {
        vector<int> count(n_bins, 0);
        for (int bin = 0; bin < (n_bins - 1); bin++)
        {
            count[bin] = sbs.db()->count_observations(source,
                                                      mjd_start + bin * step,
                                                      mjd_start + (bin + 1) * step);
            progress.update();
            std::cerr << "    \r";
            progress.status(false);
        }
        std::cerr << "\r            \n";

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

        // Set log level
        int log_level = sbsearch::INFO;
        if (args.verbose)
            log_level = sbsearch::DEBUG;

        SBSearch sbs(SBSearch::sqlite3, args.database, {args.log_file, log_level});
        Logger::info() << "SBSearch observation management tool." << std::endl;

        // Set log level
        if (args.verbose)
            Logger::get_logger().log_level(sbsearch::DEBUG);
        else
            Logger::get_logger().log_level(sbsearch::INFO);

        if (args.action == "add") // add data to database
        {
            if (args.file == "-")
            {
                cout << "\nReading observations from stdin.\n";
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
