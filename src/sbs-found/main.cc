#include <iostream>
#include <string>
#include <vector>
#include <boost/json.hpp>
#include <boost/program_options.hpp>

#include "config.h"
#include "date.h"
#include "logging.h"
#include "moving_target.h"
#include "sbsdb.h"
#include "sbsearch.h"
#include "cli.h"
#include "table.h"
#include "sbs-found/arguments.h"
#include "sbs-found/list.h"
#include "sbs-found/remove.h"
#include "sbs-found/summarize.h"

using namespace sbsearch;
using namespace sbsearch::cli;
using namespace sbsearch::table;
namespace json = boost::json;

using std::cerr;
using std::cout;
using std::string;
using std::to_string;
using std::vector;

namespace sbsearch::sbs_found
{
    template <typename DB>
    void sbs_found(int argc, char *argv[])
    {
        Arguments args = get_arguments(argc, argv);

        SBSearch<DB> sbs(args.database, {args.log_file, args.log_level()});

        Logger::info() << "SBSearch moving target found catalog tool." << std::endl;

        // Targets of interest, if any
        vector<MovingTarget> targets;
        if (args.input_file)
        {
            std::istream *is;
            std::ifstream inf;
            if (args.target == "-")
                is = &std::cin;
            else
            {
                inf.open(args.target);
                is = &inf;
            }

            for (string name; std::getline(*is, name);)
                if ((name.size() > 0) && (name[0] != '#'))
                    targets.push_back(sbsdb::get::moving_target(sbs.db(), name, !args.major_body));
        }
        else if (!args.target.empty())
            targets.push_back(sbsdb::get::moving_target(sbs.db(), args.target, !args.major_body));

        // list or summarize the observations of those targets
        if (args.action == "list")
            // when no targets are specified, get all found observations
            if (targets.size() == 0)
                list_found(sbs,
                           args.start_date.mjd(),
                           args.stop_date.mjd(),
                           args.sources,
                           args.output_filename,
                           args.output_format);
            else
                list_found(sbs,
                           targets,
                           args.start_date.mjd(),
                           args.stop_date.mjd(),
                           args.sources,
                           args.output_filename,
                           args.output_format);
        else if (args.action == "remove")
        {
            string prompt;
            if (targets.empty())
                prompt = ("Remove all found rows for all targets between " +
                          args.start_date.iso() + " and " + args.stop_date.iso() + "?");
            else
            {
                prompt = ("Remove all found rows for " + to_string(targets.size()) + " target" +
                          (targets.size() == 1 ? "" : "s") + " between " +
                          args.start_date.iso() + " and " + args.stop_date.iso() + "?");
            }

            if (args.force || confirm(prompt))
            {
                if (targets.empty())
                    remove_found(sbs, args.start_date.mjd(), args.stop_date.mjd());
                else
                    remove_found(sbs, targets, args.start_date.mjd(), args.stop_date.mjd());
            }
        }
        else if (args.action == "summarize")
        {
            // when no targets are specified, summarize all targets
            if (targets.empty())
                targets = sbsdb::get::all_moving_targets(sbs.db());

            summarize_found(sbs,
                            targets,
                            args.start_date.mjd(),
                            args.stop_date.mjd(),
                            args.sources,
                            args.output_filename,
                            args.output_format);
        }
    }
}

int main(int argc, char *argv[])
{
    try
    {
        // get basic CLI stuff first to tell us which flavor of database to use
        string database = sbs_found::get_arguments(argc, argv).database;
        if (database.substr(0, 8) == "postgres")
            sbs_found::sbs_found<sbsdb::Postgresql>(argc, argv);
        else
            throw DatabaseError("Cannot determine database type from string " + database);
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
