#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <curl/curl.h>

#include "config.h"
#include "cli.h"
#include "json.h"
#include "logging.h"
#include "sbsearch.h"
#include "sbsdb.h"
#include "sbs-query/arguments.h"
#include "sbs-query/fixed_target.h"
#include "sbs-query/moving_target.h"
#include "util/string.h"

using namespace sbsearch;
using namespace sbsearch::cli;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace sbsearch::sbs_query
{
    // make a list of targets from the CLI
    template <typename DB>
    vector<string> get_targets(const Arguments &args, SBSearch<DB> &sbs)
    {
        vector<string> targets;

        if (args.input_file)
        {
            std::ifstream input(args.target);
            for (string line; std::getline(input, line);)
                if ((line.size() > 0) && (line[0] != '#'))
                    targets.emplace_back(line);
        }
        else if (args.all_moving_targets)
        {
            for (const MovingTarget &moving_target : sbsdb::get::all_moving_targets(sbs.db()))
                targets.emplace_back(moving_target.designation());
        }
        else
            targets.emplace_back(args.target);

        return targets;
    }

    template <typename DB>
    void sbs_query(int argc, char *argv[])
    {
        Arguments args = get_arguments(argc, argv);

        if ((args.output_format == OutputFormat::AUTO) && !args.output_file.empty())
            args.output_format = get_output_format(args.output_file);

        SBSearch<DB> sbs(args.database, {args.log_file, args.log_level()});
        message::info("SBSearch target query tool " SBSEARCH_VERSION "\n");

        // setup moving/fixed target name array
        vector<string> targets = get_targets(args, sbs);

        // Extra console messages go to...
        std::ostream *console = get_console(args.quiet);

        // Set up output stream: file or stdout
        std::ostream *os;
        std::ofstream outf;
        if (args.output_file.empty())
            os = &cout;
        else
        {
            outf.open(args.output_file);
            os = &outf;
        }

        // fixed target search
        if (args.fixed_target)
        {
            Observations observations = fixed_target::query(targets, args, sbs, console);

            // output
            observations.format.date = args.date_format;
            observations.format.show_fov = args.show_fov;
            if (args.output_format == JSON)
                *os << json::value_from(observations);
            else
                *os << observations;
        }
        else
        {
            Founds founds = moving_target::query(targets, args, sbs, console);
            founds.ephemeris_format.date = args.date_format;
            founds.observation_format.date = args.date_format;
            founds.observation_format.show_fov = args.show_fov;

            // output
            if (args.save)
                cout << founds.size() << " results saved to the database." << endl;
            else
            {
                cout << "\n";
                if (args.output_format == JSON)
                    *os << json::value_from(founds);
                else
                    *os << founds;
            }
        }

        *os << "\n";

        if (!args.info_file.empty())
        {
            std::ofstream outf(args.info_file);
            outf << sbs.query_info();
        }
    }
}

int main(int argc, char *argv[])
{
    curl_global_init(CURL_GLOBAL_ALL);

    try
    {
        // identify which flavor of database to use
        string database = sbsearch::sbs_query::get_arguments(argc, argv).database;
        if (database.substr(0, 8) == "postgres")
            sbsearch::sbs_query::sbs_query<sbsdb::Postgresql>(argc, argv);
        else
            throw DatabaseError("Cannot determine database type from string " + database);
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }

    return 0;
}