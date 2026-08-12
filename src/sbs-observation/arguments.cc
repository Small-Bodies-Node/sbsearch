#include <iostream>
#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "./arguments.h"

using namespace sbsearch;
using namespace sbsearch::cli;
using std::cerr;
using std::cout;
using std::endl;
using std::optional;
using std::string;
using std::vector;

namespace sbsearch::sbs_observation
{
    Arguments get_arguments(int argc, char *argv[])
    {
        using namespace boost::program_options;

        Arguments args;

        positional_options_description positional;
        positional.add("action", 1).add("file", 2);

        options_description hidden("Hidden options");
        hidden.add_options()(
            "action", value<string>(&args.action), "action")(
            "file", value<string>(&args.file), "read data from this file");

        options_description add_options("Options for add/update actions");
        add_options.add_options()(
            "format-help", "display help on input file format and exit")(
            "format", value<OutputFormat>(&args.file_format)->default_value(OutputFormat::AUTO), "input file format: json or csv")(
            ",i", bool_switch(&args.drop_indices), "drop observations indices before adding add, re-build indices when done")(
            "batch-size,b", value<int>(&args.batch_size)->default_value(10000), "expect up to <n> observations per JSON object, or parse <n> observations at a time from CSV files")(
            ",n", bool_switch(&args.noop), "no-op mode, parse the file, but do not add to the database");

        options_description source_options("Options for data sources");
        source_options.add_options()(
            "source,s", value<vector<string>>(&args.sources), "only summarize this source, may be specified multiple times")(
            "start", value<optional<Date>>(&args.start_date),
            "start date for summarize [YYYY-MM-DD]")(
            "stop,end", value<optional<Date>>(&args.stop_date),
            "stop date for summarize [YYYY-MM-DD]");

        options_description general = get_common_options((CommonArguments *)&args);

        options_description visible("");
        visible.add(add_options).add(source_options).add(general);

        options_description all("");
        all.add(visible).add(hidden);

        options_description add_update_action("");
        add_update_action.add(add_options).add(general);

        options_description summarize_action("");
        summarize_action.add(source_options).add(general);

        variables_map vm;
        boost::program_options::store(command_line_parser(argc, argv).options(all).positional(positional).run(), vm);
        boost::program_options::notify(vm);

        if (vm.count("version"))
        {
            cout << "SBSearch " << SBSEARCH_VERSION << "\n";
            exit(0);
        }

        if (vm.count("format-help"))
        {
            cout << R"(
The CSV file format is:

  source,observatory,product_id,mjd_start,mjd_stop,fov,observation_id,meta
  "Big Survey Project","I41","unique product ID",60000.00,60000.01,"0:0, 1:0, 1:1, 0:1",1,"{\"maglim\": 25.0}"

The JSON file format is a list of objects:

  [
    {
      "source": "Big Survey Project",
      "observatory": "I41",
      "product_id": "unique product ID",
      "mjd_start": 60000.00,
      "mjd_stop": 60000.01,
      "fov": "0:0, 1:0, 1:1, 0:1",
      "observation_id": 1,
      "meta": {...}
    },
    {
      ... the next observation
    },
    ... and more
  ]
  [ ... additional arrays as needed ]

Notes:
* "fov" is a string of comma-separated RA:Dec pairs in units of degrees.
* "observation_id" is optional, but if included and it matches a record
  in the database then the database entry is updated.
* "meta" is optional.
  - Empty strings in CSV are treated as null.
  - JSON objects are converted to strings.
* Multiple JSON arrays of observations may be included in the file.  This
  helps with memory conservation when >>10000 observations are being added.
* The CSV column order is arbitrary.
* The CSV reader ignores empty lines and lines starting with #.

)";
            exit(0);
        }

        if (vm.count("help") || !vm.count("action"))
        {
            // help for a specific action?
            if (args.action == "add")
            {
                cout << "Usage: sbs-observation add <file> [options...]\n"
                     << "Add observations to the database.\n\n"
                     << "<file> contains JSON- or CSV-formatted data\n"
                     << add_update_action << "\n";
            }
            else if (args.action == "update")
            {
                cout << "Usage: sbs-observation update <file> [options...]\n"
                     << "Add observations to the database, updating as needed.\n\n"
                     << "<file> contains JSON- or CSV-formatted data\n"
                     << add_update_action << "\n";
            }
            else if (args.action == "summarize")
            {
                cout << "Usage: sbs-observation summarize [options...]\n"
                     << "Summarize the observation database.\n\n"
                     << summarize_action << "\n";
            }
            else
            {
                cout << "Usage: sbs-observation <action> [options...]\n\n"
                     << "Manage sbsearch observations.\n\n"
                     << "<action> is one of {add, update, summarize}\n"
                     << visible << "\n";
            }

            if (!vm.count("action") && !vm.count("help"))
                cout << "\naction is a required argument\n";

            exit(0);
        }

        validate_common_options(vm);
        action_is(args.action, {"add", "update", "summarize"});
        action_dependency(vm, "add", "file");

        if ((args.file == "-") && (args.file_format == AUTO))
            throw std::runtime_error("When reading from stdin, use --format.");

        return args;
    }
}