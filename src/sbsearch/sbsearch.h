#ifndef SBSEARCH_H_
#define SBSEARCH_H_

#include <optional>
#include <string>
#include <vector>
#include <s2/s2cap.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>

#include "../ephemeris.h"
#include "../exceptions.h"
#include "../found.h"
#include "../logging.h"
#include "../indexer.h"
#include "../intersection.h"
#include "../observation.h"
#include "../query_info.h"
#include "../sbsdb/find.h"

using std::optional;
using std::string;
using std::vector;

namespace sbsearch
{
    template <typename SBSDB>
    class SBSearch
    {
    public:
        // Options:
        //   - log file name
        //   - logging level
        //   - create database if it does not exist?
        //   - save query info; retrieve it later with info()
        struct Options
        {
            std::string log_file = "/dev/null";
            int log_level = sbsearch::LogLevel::INFO;
            bool create = false;
        };

        // options for find_* methods
        struct FindOptions
        {
            // Search between mjd_start and mjd_stop.
            double mjd_start = 0;
            double mjd_stop = 100000;

            // Search this data source, or all sources if not defined.
            optional<string> source;

            // Flag to account for parallax.
            bool parallax = false;

            // New properties

            // Flag to save found ephemeris results to the database.
            bool save = false;

            // Maximum number of query cells to generate.
            uint max_spatial_query_cells = 8;

            // Expand the query to cover this distance around the region.
            double padding = 0;

            // Split ephemerides into segments of this length (deg) and time period (day).
            double arc_length = 4;
            double time_period = 30;

            // Type of intersections that result in a match for fixed region queries.
            IntersectionType intersection_type = IntersectsArea;

            // return approximate results?
            bool approximate = false;

            //  save query info; retrieve it later with info()
            bool save_info = false;

            // Validate parameters
            void validate() const
            {
                if (mjd_start > mjd_stop)
                    throw SBSException("Find start date is after stop date.");
            };

            // Convert to an FindOptions object.
            sbsdb::find::Options as_sbsearch_db_options() const
            {
                return sbsdb::find::Options{mjd_start, mjd_stop, source};
            }
        };

        // constructor
        //
        // Setting options.log_file has no effect if the Logger has already been
        // initalized.
        //
        // `uri` is used to initialize the Database object.
        SBSearch(const string &uri, const Options &options = Options());

        // database maintainence
        //
        // drop/create indices, generally used when adding many new observations
        inline void drop_observations_indices() { db_.drop_observations_indices(); };
        inline void create_observations_indices() { db_.create_observations_indices(); };

        // read-only access to indexer options
        const Indexer::Options &indexer_options() { return indexer_.options(); };

        // Re-index the terms for each observation and ephemeris, and
        // store the new indexer parameters to the database.
        void reindex_database_terms(const Indexer::Options &options);

        // database I/O

        SBSDB *db() { return &db_; }

        // Add ephemeris data to the database.
        //
        // If the ephemeris's target is not already in the database, then it
        // will be added and eph.target() updated.
        //
        // If there is ephemeris data already for this target and date range,
        // then an EphemerisError is thrown.
        void add_ephemeris(Ephemeris &eph);

        // Add indices to these observations.
        void index_observations(Observations &observations);

        // Add observations, index terms will be added as needed.  Generally
        // users will use this instead of sbsdb::add::observations().
        void add_observations(Observations &observations);

        // Update observations by `product_id`, index terms will be added as
        // needed and `observation_id` may be updated.  Generally users will use
        // this instead of sbsdb::update::observations().
        void update_observations(Observations &observations);

        // search functions

        // Search for observations by point.
        Observations find_observations(const S2Point &point, const FindOptions &options = FindOptions());

        // Search for observations by spherical cap.
        Observations find_observations(const S2Cap &cap, const FindOptions &options = FindOptions());

        // Search for observations by polygon.
        Observations find_observations(const S2Polygon &polygon, const FindOptions &options = FindOptions());

        // Search for observations by ephemeris.
        Founds find_observations(const Ephemeris &ephemeris, const FindOptions &options = FindOptions());

        // Retrieve query info saved during a find_observations call with options.save_info = true.
        const QueryInfo query_info()
        {
            return query_info_;
        };

    private:
        SBSDB db_;
        Indexer indexer_;
        S2RegionTermIndexer center_indexer_{};
        QueryInfo query_info_;
    };
}

#endif // SBSEARCH_H_
