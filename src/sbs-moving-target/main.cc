#include <iostream>
#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "config.h"
#include "date.h"
#include "logging.h"
#include "moving_target.h"
#include "sbs-moving-target/arguments.h"
#include "sbs-moving-target/add.h"
#include "sbs-moving-target/list.h"
#include "sbs-moving-target/remove.h"
#include "sbs-moving-target/summarize.h"
#include "sbsdb.h"
#include "sbsearch.h"
#include "table.h"
#include "util/string.h"

using namespace sbsearch;
using namespace sbsearch::table;
using sbsearch::SBSearch;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace sbsearch::sbs_moving_target
{
    template <typename DB>
    void sbs_moving_target(int argc, char *argv[])
    {
        Arguments args = get_arguments(argc, argv);

        SBSearch<DB> sbs(args.database, {args.log_file, args.log_level()});
        Logger::info() << "SBSearch moving target management tool." << endl;

        if (args.action == "add")
            add(args.target, args.major_body, args.alternate_names, sbs);
        if (args.action == "list")
            list(sbs);
        else if (args.action == "remove")
            remove(args.target, args.major_body, args.force_remove, sbs);
        else if (args.action == "summarize")
            summarize(args.start_date, args.stop_date, sbs);
    }
}

int main(int argc, char *argv[])
{
    try
    {
        // get basic CLI stuff first to tell us which flavor of database to use

        string database = sbs_moving_target::get_arguments(argc, argv).database;
        if (database.substr(0, 8) == "postgres")
            sbs_moving_target::sbs_moving_target<sbsdb::Postgresql>(argc, argv);
        else
            throw DatabaseError("Cannot determine database type from string " + database);
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
