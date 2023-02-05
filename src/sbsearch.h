#ifndef SBSEARCH_H_
#define SBSEARCH_H_

#include <string>
#include <vector>
#include <s2/s2latlng.h>
#include <s2/s2region.h>

#include "ephemeris.h"
#include "indexer.h"
#include "sbsdb.h"

namespace sbsearch
{
    struct Found
    {
        Observation observation;
        Ephemeris ephemeris;

        Found(Observation o, Ephemeris e) : observation(o), ephemeris(e){};
    };

    class SBSearch
    {
    public:
        enum DatabaseType
        {
            sqlite3
        };

        // constructor
        //
        // For sqlite3 databases:
        //   - `name` is the database filename, ":memory:" for an in-memory
        //     database, or "" (empty-string) for a temporary on-disk database.
        SBSearch(DatabaseType database_type, const char *name, Indexer::Options indexer_options = Indexer::Options());

        void add_observations(vector<Observation> &observations);

        // search functions
        //
        // Search by point or polygon, optionally over a time range.  Let
        // `mjd_start` or `mjd_stop` = -1 for an unbounded limit.
        vector<Observation> find_observations(const S2Point &point, double mjd_start = -1, double mjd_stop = -1);
        vector<Observation> find_observations(const S2Polygon &polygon, double mjd_start = -1, double mjd_stop = -1);

        // Search by ephemeris.
        vector<Found> find_observations(const Ephemeris &ephemeris);

    private:
        SBSearchDatabase *db_;
        Indexer indexer_;
    };
}

#endif // SBSEARCH_H_
