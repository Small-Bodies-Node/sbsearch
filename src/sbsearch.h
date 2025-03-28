#ifndef SBSEARCH_H_
#define SBSEARCH_H_

#include <string>
#include <vector>
#include <s2/s2latlng.h>
#include <s2/s2region.h>

#include "ephemeris.h"
#include "found.h"
#include "logging.h"
#include "indexer.h"
#include "observation.h"
#include "observatory.h"
#include "sbsdb/find.h"

namespace sbsearch
{
    // Intersection types
    enum IntersectionType
    {
        ContainsPoint,
        ContainsArea,
        IntersectsArea,
        ContainedByArea,
        ContainsCenter = ContainsPoint
    };

    std::istream &operator>>(std::istream &in, IntersectionType &intersection_type);

    template <typename SBSDB>
    class SBSearch
    {
    public:
        // Options:
        //   - log file name
        //   - logging level
        //   - create database if it does not exist?
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

            // Search this data source, or all sources if empty.
            string source = string();

            // Flag to account for parallax.
            bool parallax = false;

            // New properties

            // Flag to save found ephemeris results to the database.
            bool save = false;

            // Maximum number of query cells to generate.
            int max_spatial_query_cells = 8;

            // Expand the query to cover this distance around the region.
            double padding = 0;

            // Split ephemerides into segments of this length (deg) and time period (yr).
            double arc_length = 10;
            double time_period = 1;

            // Type of intersections that result in a match for fixed region queries.
            IntersectionType intersection_type = IntersectsArea;

            // return approximate results?
            bool approximate = false;

            // Convert to an FindOptions object.
            sbsdb::find::Options as_sbsearch_db_options() const
            {
                return sbsdb::find::Options{mjd_start,
                                            mjd_stop,
                                            source,
                                            parallax};
            }
        };

        // constructor
        //
        // Setting options.log_file has no effect if the Logger has already been
        // initalized.
        //
        // `uri` is used to initialize the Database object.
        SBSearch(const string &uri, const Options &options = Options()) : db_(uri)
        {
            // attempt to initialize logger
            Logger::get_logger(options.log_file).log_level(options.log_level);

            if (options.create)
                db_.setup_tables();

            indexer_ = Indexer(db_.indexer_options());
        };

        // database maintainence
        //
        // drop/create indices, generally used when adding many new observations
        inline void drop_observations_indices() { db_.drop_observations_indices(); };
        inline void create_observations_indices() { db_.create_observations_indices(); };

        // read-only access to indexer options
        const Indexer::Options &indexer_options() { return indexer_.options(); };

        // Re-index the terms for each observation and ephemeris, and
        // store the new indexer parameters to the database.
        void reindex(Indexer::Options options);

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

        // Add observations, index terms will be added as needed.  Generally
        // users will use this instead of db()->add_observations().
        void add_observations(Observations &observations);

        // search functions

        // Search for observations by point.
        Observations find_observations(const S2Point &point, const FindOptions &options = FindOptions());

        // Search for observations by polygon.
        Observations find_observations(const S2Polygon &polygon, const FindOptions &options = FindOptions());

        // Search for observations by ephemeris.
        Founds find_observations(const Ephemeris &ephemeris, const FindOptions &options = FindOptions());

        // Test for intersection between a polygon and a spherical cap.
        static bool intersects(const S2Polygon &polygon, const S2Cap &area, const IntersectionType intersection_type);

        // Test for intersection between two polygons.
        static bool intersects(const S2Polygon &polygon, const S2Polygon &area, const IntersectionType intersection_type);

    private:
        SBSDB db_;
        Indexer indexer_;
    };

    // std::istream &operator>>(std::istream &in, SBSearch::DatabaseType &db_type);
}

#endif // SBSEARCH_H_
