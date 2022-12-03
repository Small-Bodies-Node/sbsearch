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

#define UNDEFINED_OBSID -1

namespace sbsearch
{
    class Observation
    {
    public:
        // Initialize from values
        Observation(double mjd_start, double mjd_stop, const char *fov, int64 observation_id = UNDEFINED_OBSID);
        Observation(double mjd_start, double mjd_stop, S2LatLngRect fov, int64 observation_id = UNDEFINED_OBSID);
        Observation(double mjd_start, double mjd_stop, double *fov, int64 observation_id = UNDEFINED_OBSID);

        // Initialize from sbsearch sqlite3 database
        // Observation(sqlite3 *db, int64 observation_id);

        // Property access
        inline int64 observation_id() { return observation_id_; };
        inline double mjd_start() { return mjd_start_; };
        inline double mjd_stop() { return mjd_stop_; };
        inline string fov() { return string(fov_); };

        // check if observation is valid
        bool is_valid();

        // observation IDs may be updated if they are not already defined
        inline void observation_id(int64 observation_id) { observation_id_ = observation_id; };

        // Generate index terms
        vector<string> index_terms(S2RegionTermIndexer &indexer);

        // Return an S2Polygon describing this observation's field-of-view
        unique_ptr<S2Polygon> as_polygon();

    private:
        int64 observation_id_;
        double mjd_start_;
        double mjd_stop_;
        char fov_[512];

        void copy_fov(char *fov);
    };
}
#endif // OBSERVATION_H_