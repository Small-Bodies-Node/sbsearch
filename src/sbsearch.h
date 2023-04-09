#ifndef SBSEARCH_H_
#define SBSEARCH_H_

#include <string>
#include <vector>
#include <s2/s2latlng.h>
#include <s2/s2region.h>

#include "ephemeris.h"
#include "found.h"
#include "indexer.h"
#include "observation.h"
#include "observatory.h"
#include "sbsdb.h"

namespace sbsearch
{
    std::ostream &operator<<(std::ostream &os, const Found &found);
    // If show_fov is set on any observation, then the FOV is shown for all
    std::ostream &operator<<(std::ostream &os, const vector<Found> &founds);

    struct SBSearchOptions
    {
        std::string log_file = "/dev/null";
        int log_level = sbsearch::INFO;
        bool create = false;
    };

    class SBSearch
    {
    public:
        enum DatabaseType
        {
            sqlite3
        };

        using Options = SBSearchOptions;

        // constructor
        //
        // Setting log_file has no effect if the Logger has already been initalized.
        //
        // For sqlite3 databases:
        //   - `name` is the database filename, ":memory:" for an in-memory
        //     database, or "" (empty-string) for a temporary on-disk database.
        SBSearch(DatabaseType database_type, const std::string name, const Options options = Options());

        ~SBSearch() { db_->close(); }

        // database maintainence
        //
        // drop/create indices, generally used when adding many new observations
        inline void drop_observations_indices() { db_->drop_observations_indices(); };
        inline void create_observations_indices() { db_->create_observations_indices(); };

        // read-only access to indexer options
        const Indexer::Options &indexer_options() { return indexer_.options(); };

        // Re-index the terms for each observation and ephemeris, and
        // store the new indexer parameters to the database.
        void reindex(Indexer::Options options);

        // database I/O

        // Most user ops can use const access to db.
        const SBSearchDatabase *db() { return db_; }

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

        // Get observations by observation ID vector.
        Observations get_observations(const vector<int64> &observation_id);

        // search options
        typedef SBSearchDatabase::Options SearchOptions;

        // search functions

        // Search for observations by point.
        Observations find_observations(const S2Point &point, const SearchOptions &options = SearchOptions());

        // Search for observations by polygon.
        Observations find_observations(const S2Polygon &polygon, const SearchOptions &options = SearchOptions());

        // Search for observations by ephemeris.
        vector<Found> find_observations(const Ephemeris &ephemeris, const SearchOptions &options = SearchOptions());

    private:
        SBSearchDatabase *db_;
        Indexer indexer_;
    };
}

#endif // SBSEARCH_H_
