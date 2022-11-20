#ifndef OBSERVATION_H_
#define OBSERVATION_H_

#include <string>
#include <vector>
#include <sqlite3.h>
#include "s2/s2polygon.h"
#include "s2/s2latlng_rect.h"
#include <s2/s2region_term_indexer.h>

using std::string;
using std::unique_ptr;
using std::vector;

namespace sbsearch
{
    class Observation
    {
    public:
        // Initialize from values
        Observation(int64 obsid, double mjd_start, double mjd_stop, char *fov);
        Observation(int64 obsid, double mjd_start, double mjd_stop, S2LatLngRect fov);
        Observation(int64 obsid, double mjd_start, double mjd_stop, double *fov);

        // Initialize from sbsearch sqlite3 database
        Observation(sqlite3 *db, int64 obsid);

        // Property access
        inline int64 obsid() { return obsid_; };
        inline double mjd_start() { return mjd_start_; };
        inline double mjd_stop() { return mjd_stop_; };
        inline string fov() { return string(fov_); };

        // Generate index terms
        vector<string> index_terms(S2RegionTermIndexer &indexer);

        // SQL statement to add to database
        void add_sql(S2RegionTermIndexer &indexer, char *statement);

        // Add to database
        void add_to_database(sqlite3 *db, S2RegionTermIndexer &indexer);

        // Return an S2Polygon describing this observation's field-of-view
        unique_ptr<S2Polygon> as_polygon();

    private:
        int64 obsid_;
        double mjd_start_;
        double mjd_stop_;
        char fov_[512];

        void copy_fov(char *fov);
    };
}
#endif // OBSERVATION_H_