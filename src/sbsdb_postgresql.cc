#include "config.h"

#include <algorithm>
#include <ctime>
#include <iostream>
#include <stdexcept>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <pqxx/pqxx>

#include "ephemeris.h"
#include "found.h"
#include "exceptions.h"
#include "logging.h"
#include "moving_target.h"
#include "observation.h"
#include "observatory.h"
#include "sbsdb.h"
#include "sbsdb_postgresql.h"

using std::endl;
using std::set;
using std::string;
using std::vector;

namespace sbsearch
{
    SBSearchDatabasePostgreSQL::SBSearchDatabasePostgreSQL(const string uri) : connection_(uri)
    {
        if (!connection_.is_open())
            throw "Failed to open database connection.";
        // try
        // {
        //     // Start a transaction.  In libpqxx, you always work in one.
        //     pqxx::work tx(connection);

        //     // work::exec1() executes a query returning a single row of data.
        //     // We'll just ask the database to return the number 1 to us.
        //     pqxx::row r = tx.exec1("SELECT 1");

        //     // Commit your transaction.  If an exception occurred before this
        //     // point, execution will have left the block, and the transaction will
        //     // have been destroyed along the way.  In that case, the failed
        //     // transaction would implicitly abort instead of getting to this point.
        //     tx.commit();

        //     // Look at the first and only field in the row, parse it as an integer,
        //     // and print it.
        //     //
        //     // "r[0]" returns the first field, which has an "as<...>()" member
        //     // function template to convert its contents from their string format
        //     // to a type of your choice.
        //     std::cout << r[0].as<int>() << std::endl;
        // }
        // catch (std::exception const &e)
        // {
        //     std::cerr << e.what() << std::endl;
        //     return;
        // }
    }

    void SBSearchDatabasePostgreSQL::close()
    {
        if (connection_.is_open())
        {
            connection_.close();
            Logger::info() << "Closed database." << endl;
        }
    };

    void SBSearchDatabasePostgreSQL::setup_tables() {

    };

    void SBSearchDatabasePostgreSQL::drop_observations_indices() {};
    void SBSearchDatabasePostgreSQL::create_observations_indices() {};
    void SBSearchDatabasePostgreSQL::execute_sql(const char *statement) const {};

    // get single value results from a SQL statement
    double *SBSearchDatabasePostgreSQL::get_double(const char *statement) const {};
    int *SBSearchDatabasePostgreSQL::get_int(const char *statement) const {};
    int64 *SBSearchDatabasePostgreSQL::get_int64(const char *statement) const {};
    string *SBSearchDatabasePostgreSQL::get_string(const char *statement) const {};

    void SBSearchDatabasePostgreSQL::indexer_options(Indexer::Options options) {};

    std::pair<double *, double *> SBSearchDatabasePostgreSQL::observation_date_range(const string &source) const {};

    void SBSearchDatabasePostgreSQL::add_moving_target(MovingTarget &target) const {};
    void SBSearchDatabasePostgreSQL::remove_moving_target(const MovingTarget &target) const {};
    void SBSearchDatabasePostgreSQL::update_moving_target(const MovingTarget &target) const {};
    MovingTarget SBSearchDatabasePostgreSQL::get_moving_target(const int moving_target_id) const {};
    MovingTarget SBSearchDatabasePostgreSQL::get_moving_target(const string &name, const bool small_body) const {};
    vector<MovingTarget> SBSearchDatabasePostgreSQL::get_all_moving_targets() const {};

    void SBSearchDatabasePostgreSQL::add_observatory(const string &name, const Observatory &observatory) const {};
    const Observatory SBSearchDatabasePostgreSQL::get_observatory(const string &name) const {};
    const Observatories SBSearchDatabasePostgreSQL::get_observatories() const {};
    void SBSearchDatabasePostgreSQL::remove_observatory(const string &name) const {};
    const vector<string> SBSearchDatabasePostgreSQL::get_sources() const {};

    void SBSearchDatabasePostgreSQL::add_ephemeris(Ephemeris &eph) const {};
    Ephemeris SBSearchDatabasePostgreSQL::get_ephemeris(const MovingTarget target, double mjd_start, double mjd_stop) const {};
    int SBSearchDatabasePostgreSQL::remove_ephemeris(const MovingTarget target, double mjd_start, double mjd_stop) const {};
    std::pair<double *, double *> SBSearchDatabasePostgreSQL::ephemeris_date_range() const {};

    void SBSearchDatabasePostgreSQL::add_observation(Observation &observation) const {};
    Observation SBSearchDatabasePostgreSQL::get_observation(const int64 observation_id) const {};
    void SBSearchDatabasePostgreSQL::remove_observations(const double mjd_start, const double mjd_stop) const {};
    void SBSearchDatabasePostgreSQL::remove_observations(const string &source, const double mjd_start, const double mjd_stop) const {};

    // Count number of observations within an interval.
    int64 SBSearchDatabasePostgreSQL::count_observations(const double mjd_start, const double mjd_stop) const {};

    // Count number of observations for a source within an interval, if
    // source is an empty string, then count all sources.
    int64 SBSearchDatabasePostgreSQL::count_observations(const string &source, const double mjd_start, const double mjd_stop) const {};

    Observations SBSearchDatabasePostgreSQL::find_observations(const double mjd_start, const double mjd_stop, const int64 limit, const int64 offset) const {};
    Observations SBSearchDatabasePostgreSQL::find_observations(const string &source, const double mjd_start, double mjd_stop, const int64 limit, const int64 offset) const {};
    Observations SBSearchDatabasePostgreSQL::find_observations(vector<string> query_terms, const Options &options) const {};

    void SBSearchDatabasePostgreSQL::add_found(const Found &found) const {};
    Founds SBSearchDatabasePostgreSQL::get_found(const Observation &observation) const {};
    Founds SBSearchDatabasePostgreSQL::get_found(const MovingTarget &target) const {};
    void SBSearchDatabasePostgreSQL::remove_found(const Found &found) const {};

}