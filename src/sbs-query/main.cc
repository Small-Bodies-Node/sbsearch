#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <curl/curl.h>

#include "config.h"
#include "logging.h"
#include "sbsearch.h"
#include "sbsdb.h"
#include "cli.h"
#include "sbs-query/arguments.h"
#include "sbs-query/fixed_target.h"
#include "sbs-query/moving_target.h"

using namespace sbsearch;
using namespace sbsearch::cli;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace sbsearch::sbs_query
{
    template <typename DB>
    void sbs_query(int argc, char *argv[])
    {
        Arguments args = get_arguments(argc, argv);

        if ((args.output_format == OutputFormat::AUTO) && !args.output_file.empty())
            args.output_format = get_output_format(args.output_file);

        SBSearch<DB> sbs(args.database, {args.log_file, args.log_level()});
        message("SBSearch target query tool.\n");

        // setup moving/fixed target name array
        vector<string> targets;
        if (args.input_file)
        {
            std::ifstream input(args.target);
            for (string line; std::getline(input, line);)
                if ((line.size() > 0) && (line[0] != '#'))
                    targets.push_back(line);
        }
        else if (args.all_moving_targets)
        {
            for (const MovingTarget &moving_target : sbsdb::get::all_moving_targets(sbs.db()))
                targets.push_back(moving_target.designation());
        }
        else
            targets.push_back(args.target);

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
            Observations observations;
            for (string target : targets)
            {
                Observations new_observations = query_fixed_target(args, target, sbs);
                observations.append(new_observations);
            }

            // output
            observations.format.show_fov = args.show_fov;
            if (args.output_format == JSON)
            {
                json::array array;
                for (Observation obs : observations)
                    array.emplace_back(obs.as_json());

                *os << array;
            }
            else
            {
                *os << observations;
            }
        }
        else
        // moving target search
        {
            // this works for searches for single targets by ephemeris file or orbit
            // because get_arguments() prevents multiple targets when eph_file or
            // orbit_file are specified
            Founds founds;
            for (string target : targets)
                founds.append(query_moving_target(args, target, sbs));

            // output, but only when not saving to the database
            if (!args.save)
            {
                cout << "\n";
                founds.ephemeris_format.date = args.date_format;
                founds.observation_format.show_fov = args.show_fov;
                if (args.output_format == JSON)
                    *os << founds.as_json();
                else
                    *os << founds;
            }
            else
            {
                cout << founds.size() << " results saved to the database." << endl;
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