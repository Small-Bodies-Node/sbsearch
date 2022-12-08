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
        Observation(double mjd_start, double mjd_stop, string fov, string terms = "", int64 observation_id = UNDEFINED_OBSID);
        Observation(double mjd_start, double mjd_stop, vector<S2LatLng> vertices, string terms = "", int64 observation_id = UNDEFINED_OBSID);

        // Property access
        inline int64 observation_id() { return observation_id_; };
        inline double mjd_start() { return mjd_start_; };
        inline double mjd_stop() { return mjd_stop_; };
        inline string fov() { return string(fov_); };
        inline string terms() { return string(terms_); };

        // check if observation is valid
        bool is_valid();

        // test if observation has the same FOV as another
        bool is_same_fov(Observation &other);

        // test if observation is equal to another by comparing
        // - observation_id
        // - mjd_start
        // - mjd_stop
        // - fov
        bool is_equal(Observation &other);

        // observation IDs may be updated if they are not already defined
        inline void observation_id(int64 observation_id)
        {
            if (observation_id_ != UNDEFINED_OBSID)
                throw std::runtime_error("Observation ID already defined.");
            else
                observation_id_ = observation_id;
        };

        // Generate index terms from indexer, optionally update the terms property
        vector<string> index_terms(S2RegionTermIndexer &indexer, bool update = true);

        // Generate query terms from indexer
        vector<string> query_terms(S2RegionTermIndexer &indexer);

        // Return an S2Polygon describing this observation's field-of-view
        unique_ptr<S2Polygon> as_polygon();

    private:
        int64 observation_id_;
        double mjd_start_;
        double mjd_stop_;
        string fov_;
        string terms_;

        enum TermStyle
        {
            index,
            query
        };
        vector<string> generate_terms(TermStyle style, S2RegionTermIndexer &indexer);
    };
}
#endif // OBSERVATION_H_