#include <iostream>
#include <string>
#include <vector>
#include <curl/curl.h>

#include "config.h"
#include "ephemeris.h"
#include "logging.h"
#include "sbsearch.h"
#include "sbsdb/get.h"
#include "sbsdb/postgresql.h"
#include "sbs-ephemeris/add.h"
#include "sbs-ephemeris/arguments.h"
#include "sbs-ephemeris/get.h"
#include "sbs-ephemeris/list.h"
#include "sbs-ephemeris/remove.h"

using namespace sbsearch;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace sbsearch::sbs_ephemeris
{
    template <typename DB>
    void sbs_ephemeris(int argc, char *argv[])
    {
        Arguments args = get_arguments(argc, argv);

        SBSearch<DB> sbs(args.database, {args.log_file, args.log_level()});

        Logger::info() << "SBSearch ephemeris management tool." << endl;

        vector<string> targets;
        if (args.input_file)
        {
            std::ifstream input(args.target);
            for (string line; std::getline(input, line);)
                if ((line.size() > 0) && (line[0] != '#'))
                    targets.push_back(line);
        }
        else
            targets.push_back(args.target);

        for (string target : targets)
        {
            args.target = target;
            if (args.action == "add") // add data to database
                add(args, sbs);
            if (args.action == "get") // get ephemeris from Horizons
                get(args, sbs);
            else if (args.action == "list") // list ephemeris data
                list(args, sbs);
            else if (args.action == "remove") // remove data from database
                remove(args, sbs);
        }
    }
}

int main(int argc, char *argv[])
{
    curl_global_init(CURL_GLOBAL_ALL);

    try
    {
        // get basic CLI stuff first to tell us which flavor of database to use
        string database = sbsearch::sbs_ephemeris::get_arguments(argc, argv).database;
        if (database.substr(0, 8) == "postgres")
            sbsearch::sbs_ephemeris::sbs_ephemeris<sbsdb::Postgresql>(argc, argv);
        else
            throw DatabaseError("Cannot determine database type from string " + database);
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
