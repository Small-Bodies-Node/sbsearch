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
        Observation() = delete;

        // Initialize from values
        Observation(double mjd_start, double mjd_stop, const char *fov, string terms = "", int64 observation_id = UNDEFINED_OBSID);
        Observation(double mjd_start, double mjd_stop, double *lat, double *lng, string terms = "", int64 observation_id = UNDEFINED_OBSID)
            : Observation(mjd_start, mjd_stop, "", terms, observation_id)
        {
            sprintf(fov_, "%f:%f, %f:%f, %f:%f, %f:%f",
                    lat[0], lng[0],
                    lat[1], lng[1],
                    lat[2], lng[2],
                    lat[3], lng[3]);
        }
        Observation(double mjd_start, double mjd_stop, S2LatLngRect fov, string terms = "", int64 observation_id = UNDEFINED_OBSID)
            : Observation(mjd_start, mjd_stop, "", terms, observation_id)
        {
            sprintf(fov_, "%f:%f, %f:%f, %f:%f, %f:%f",
                    fov.lat_lo().degrees(), fov.lng_lo().degrees(),
                    fov.lat_lo().degrees(), fov.lng_hi().degrees(),
                    fov.lat_hi().degrees(), fov.lng_hi().degrees(),
                    fov.lat_hi().degrees(), fov.lng_lo().degrees());
        }

        // Initialize from sbsearch sqlite3 database
        // Observation(sqlite3 *db, int64 observation_id);

        // Property access
        inline int64 observation_id() { return observation_id_; };
        inline double mjd_start() { return mjd_start_; };
        inline double mjd_stop() { return mjd_stop_; };
        inline string fov() { return string(fov_); };
        inline string terms() { return string(terms_); };

        // check if observation is valid
        bool is_valid();

        // observation IDs may be updated if they are not already defined
        inline void observation_id(int64 observation_id)
        {
            if (observation_id_ == UNDEFINED_OBSID)
                throw std::runtime_error("Observation ID already defined.");
            else
                observation_id_ = observation_id;
        };

        // Generate index terms from indexer
        void terms(S2RegionTermIndexer &indexer);

        // Return an S2Polygon describing this observation's field-of-view
        unique_ptr<S2Polygon> as_polygon();

    private:
        int64 observation_id_;
        double mjd_start_;
        double mjd_stop_;
        char fov_[512];
        string terms_;

        void copy_fov(char *fov);
    };
}
#endif // OBSERVATION_H_