#ifndef OBSERVATION_H_
#define OBSERVATION_H_

#include <functional>
#include <string>
#include <vector>
#include <s2/s2polygon.h>

#include "util.h"

using sbsearch::format_vertices;
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
        Observation(string source, string product_id, double mjd_start, double mjd_stop, string fov, string terms = "", int64 observation_id = UNDEFINED_OBSID);
        Observation(string source, string product_id, double mjd_start, double mjd_stop, vector<S2LatLng> vertices, string terms = "", int64 observation_id = UNDEFINED_OBSID)
            : Observation(source, product_id, mjd_start, mjd_stop, format_vertices(vertices), terms, observation_id){};

        // Property getters
        inline string source() const { return source_; };
        inline string product_id() const { return product_id_; };
        inline int64 observation_id() const { return observation_id_; };
        inline double mjd_start() const { return mjd_start_; };
        inline double mjd_stop() const { return mjd_stop_; };
        inline string fov() const { return string(fov_); };
        inline string terms() const { return string(terms_); };

        // Property setters
        inline void source(const string new_source) { source_ = string(new_source); };
        inline void product_id(const string new_product_id) { product_id_ = string(new_product_id); };
        void observation_id(int64 new_observation_id);
        inline void mjd_start(double new_mjd_start) { mjd_start_ = new_mjd_start; };
        inline void mjd_stop(double new_mjd_stop) { mjd_stop_ = new_mjd_stop; };
        inline void fov(string new_fov) { fov_ = string(new_fov); };
        void terms(const vector<string> new_terms);
        void terms(const string new_terms);

        // check if observation is valid
        bool is_valid() const;

        // test if observation has the same FOV as another
        bool is_same_fov(const Observation &other) const;

        // test if observation is equal to another by comparing
        // - source
        // - product_id
        // - mjd_start
        // - mjd_stop
        // - fov
        bool is_equal(const Observation &other) const;
        bool operator==(const Observation &other) const;

        S2Polygon as_polygon() const;

    private:
        string source_;
        string product_id_;
        int64 observation_id_;
        double mjd_start_;
        double mjd_stop_;
        string fov_;
        string terms_;
    };

}

// custom specialization of std::hash for unordered_set<Observation>
template <>
struct std::hash<sbsearch::Observation>
{
    std::size_t operator()(sbsearch::Observation const &observation) const noexcept
    {
        char s[256];
        sprintf(s, "%lld:%s:%lf:%lf", observation.observation_id(), observation.fov().c_str(),
                observation.mjd_start(), observation.mjd_stop());
        return std::hash<std::string>{}(string(s));
    }
};
#endif // OBSERVATION_H_