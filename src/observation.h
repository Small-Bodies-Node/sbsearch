#ifndef OBSERVATION_H_
#define OBSERVATION_H_

#include <string>
#include <vector>
#include <s2/s2polygon.h>

using std::string;
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

        // Property getters
        inline int64 observation_id() const { return observation_id_; };
        inline double mjd_start() const { return mjd_start_; };
        inline double mjd_stop() const { return mjd_stop_; };
        inline string fov() const { return string(fov_); };
        inline string terms() const { return string(terms_); };

        // Property setters
        void observation_id(int64 new_observation_id);
        inline void mjd_start(double new_mjd_start) { mjd_start_ = new_mjd_start; };
        inline void mjd_stop(double new_mjd_stop) { mjd_stop_ = new_mjd_stop; };
        inline void fov(string new_fov) { fov_ = string(new_fov); };
        void terms(const vector<string> new_terms);
        void terms(const string new_terms);

        // check if observation is valid
        bool is_valid() const;

        // test if observation has the same FOV as another
        bool is_same_fov(Observation &other) const;

        // test if observation is equal to another by comparing
        // - observation_id
        // - mjd_start
        // - mjd_stop
        // - fov
        bool is_equal(Observation &other) const;

        S2Polygon as_polygon() const;

    private:
        int64 observation_id_;
        double mjd_start_;
        double mjd_stop_;
        string fov_;
        string terms_;

        // enum TermStyle
        // {
        //     index,
        //     query
        // };
        // vector<string> generate_terms(TermStyle style, S2RegionTermIndexer &indexer);
    };
}
#endif // OBSERVATION_H_