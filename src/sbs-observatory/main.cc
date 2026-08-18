#include <iostream>
#include <string>
#include <boost/program_options.hpp>

#include "config.h"
#include "logging.h"
#include "observatory.h"
#include "sbsearch.h"
#include "sbsdb/add.h"
#include "sbsdb/get.h"
#include "sbsdb/postgresql.h"
#include "sbsdb/remove.h"
#include "sbs-observatory/arguments.h"
#include "sbs-observatory/list.h"

using namespace sbsearch;
using std::cerr;
using std::cout;
using std::string;

namespace sbsearch::sbs_observatory
{
    template <typename DB>
    void sbs_observatory(int argc, char *argv[])
    {
        Arguments args = get_arguments(argc, argv);

        SBSearch<DB> sbs(args.database, {args.log_file, args.log_level()});
        Logger::info() << "SBSearch observatory management tool." << std::endl;

        if (args.action == "add") // add data to database
            sbsdb::add::observatory(sbs.db(), args.observatory);
        else if (args.action == "list")
            list(args.output_filename, sbs);
        else if (args.action == "remove")
            sbsdb::remove::observatory(sbs.db(), args.observatory.name);
    }
}

int main(int argc, char *argv[])
{
    try
    {
        // get flavor of database to use
        string database = sbs_observatory::get_arguments(argc, argv).database;
        if (database.substr(0, 8) == "postgres")
            sbs_observatory::sbs_observatory<sbsdb::Postgresql>(argc, argv);
        else
            throw DatabaseError("Cannot determine database type from string " + database);
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
