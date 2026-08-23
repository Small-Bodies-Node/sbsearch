#ifndef SBSEARCH_H_
#define SBSEARCH_H_

#include <optional>
#include <string>
#include <vector>
#include <s2/s2cap.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>

#include "config.h"
#include "exceptions.h"
#include "found.h"
#include "logging.h"
#include "indexer.h"
#include "intersection.h"
#include "observation.h"
#include "observatory.h"
#include "query_info.h"
#include "queue.h"
#include "ephemeris/ephemeris.h"
#include "sbsdb/search.h"

using sbsearch::ephemeris::Ephemeris;
using std::optional;
using std::string;
using std::vector;

namespace sbsearch
{
    // options for search_* methods
    struct SearchOptions
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

        // Use ephemeris uncertainty for searches.
        bool use_ephemeris_uncertainty = false;

        // Split ephemerides into segments of this length (deg) and time period (day).
        double arc_length = 4;
        double time_period = 30;

        // Type of intersections that result in a match for fixed region queries.
        IntersectionType intersection_type = IntersectsArea;

        // return approximate results?
        bool approximate = false;

        // save query info; retrieve it later with info()
        bool save_info = false;

        // console quiet mode?
        bool quiet = false;

        // console verbose mode?
        bool verbose = false;

        // debug logging?
        bool debug = false;

        // Validate parameters
        void validate() const
        {
            if (mjd_start > mjd_stop)
                throw SBSException("Find start date is after stop date.");
        };

        // Convert to a SearchOptions object.
        sbsdb::search::Options as_sbsearch_db_options() const
        {
            return sbsdb::search::Options{mjd_start, mjd_stop, source};
        }
    };

    template <class SBSDB>
    class SBSearch
    {
    public:
        // Options:
        //   - log file name
        //   - logging level
        //   - create database if it does not exist?
        //   - number of threads for intersection testing and observation
        //     re-indexing
        struct Options
        {
            std::string log_file = "/dev/null";
            int log_level = sbsearch::LogLevel::INFO;
            bool create = false;
            unsigned char threads = 2;
        };

        // constructor
        //
        // Setting options.log_file has no effect if the Logger has already been
        // initalized.
        //
        // `uri` is used to initialize the Database object.
        SBSearch(string uri, const Options &options = Options());

        // number of CPU threads to use for re-indexing and moving target
        // intersection testing
        inline unsigned char threads() { return threads_; }

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

        // Add ephemeris data to the database.·
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
        Observations search_observations(const S2Point &point, const SearchOptions &options = SearchOptions());

        // Search for observations by spherical cap.
        Observations search_observations(const S2Cap &cap, const SearchOptions &options = SearchOptions());

        // Search for observations by polygon.
        Observations search_observations(const S2Polygon &polygon, const SearchOptions &options = SearchOptions());

        // Search for observations by ephemeris.
        Founds search_observations(const Ephemeris &ephemeris, const SearchOptions &options = SearchOptions());

        // Retrieve query info saved during a search_observations call with options.save_info = true.
        const QueryInfo &query_info() { return query_info_; };

    private:
        SBSDB db_;
        Indexer indexer_;
        S2RegionTermIndexer center_indexer_{};
        QueryInfo query_info_;
        unsigned char threads_;
    };

    struct IndexFovTaskResult
    {
        vector<int64_t> observation_id;
        vector<Observation::Terms> terms;
    };

    IndexFovTaskResult index_fov_task_(Queue<std::pair<int64_t, string>> &queue,
                                       const Indexer::Options &options);
}

#endif // SBSEARCH_H_
