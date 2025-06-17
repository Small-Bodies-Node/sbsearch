#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <boost/json.hpp>
#include <boost/system.hpp>

#include "add.h"
#include "arguments.h"
#include "summary.h"
#include "../exceptions.h"
#include "../logging.h"
#include "../sbsearch.h"
#include "../sbsdb/postgresql.h"

using namespace sbsearch::sbs_observation;
using sbsearch::sbsdb::Postgresql;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace sbsearch::sbs_observation
{
    template <typename DB>
    void sbs_observation(int argc, char *argv[])
    {
        Arguments args = get_arguments(argc, argv);

        // Set log level
        int log_level = INFO;
        if (args.verbose)
            log_level = DEBUG;

        SBSearch<DB> sbs(args.database, {args.log_file, log_level});
        Logger::info() << "SBSearch observation management tool." << endl;

        // Set log level
        if (args.verbose)
            Logger::get_logger().log_level(sbsearch::DEBUG);
        else
            Logger::get_logger().log_level(sbsearch::INFO);

        if (args.action == "add") // add data to database
            add(args, sbs);
        else if (args.action == "summary")
            summary<Postgresql>(args, sbs);
    }
}

int main(int argc, char *argv[])
{
    try
    {
        // get basic CLI stuff first to tell us which flavor of database to use
        string database = get_arguments(argc, argv).database;
        if (database.substr(0, 8) == "postgres")
            sbs_observation<Postgresql>(argc, argv);
        else
            throw sbsearch::DatabaseError("Cannot determine database type from string " + database);
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
