#ifndef OBSERVATION_H_
#define OBSERVATION_H_

#include <functional>
#include <optional>
#include <ostream>
#include <string>
#include <vector>
#include <boost/json.hpp>
#include <s2/s2polygon.h>

#include "util/polygon.h"
#include "util/string.h"

using std::optional;
using std::string;
using std::string_view;
using std::vector;
namespace json = boost::json;

namespace sbsearch
{
    class Observation
    {
    public:
        Observation() {};

        // Initialize from values
        Observation(string source,
                    string observatory,
                    string product_id,
                    double mjd_start,
                    double mjd_stop,
                    string fov,
                    vector<string> terms = {},
                    optional<int64_t> observation_id = {},
                    optional<string> center = {},
                    optional<string> meta = {});

        // Copy constructor
        Observation(const Observation &other)
            : source_(other.source_),
              observatory_(other.observatory_),
              product_id_(other.product_id_),
              mjd_start_(other.mjd_start_),
              mjd_stop_(other.mjd_stop_),
              fov_(other.fov_),
              terms_(other.terms_),
              observation_id_(other.observation_id_),
              center_(other.center_),
              meta_(other.meta_)
        {
        }

        // Move constructor
        Observation(Observation &&other)
            : source_(std::move(other.source_)),
              observatory_(std::move(other.observatory_)),
              product_id_(std::move(other.product_id_)),
              mjd_start_(std::move(other.mjd_start_)),
              mjd_stop_(std::move(other.mjd_stop_)),
              fov_(std::move(other.fov_)),
              terms_(std::move(other.terms_)),
              observation_id_(std::move(other.observation_id_)),
              center_(std::move(other.center_)),
              meta_(std::move(other.meta_))
        {
            other.source_ = "";
            other.observatory_ = "";
            other.product_id_ = "";
            other.mjd_start_ = 0;
            other.mjd_stop_ = 0;
            other.fov_ = "";
            other.terms_.clear();
            other.observation_id_ = std::nullopt;
            other.center_ = std::nullopt;
            other.meta_ = std::nullopt;
        }

        // Copy assignment
        Observation &operator=(const Observation &other)
        {
            if (this != &other)
            {
                this->source_ = other.source_;
                this->observatory_ = other.observatory_;
                this->product_id_ = other.product_id_;
                this->mjd_start_ = other.mjd_start_;
                this->mjd_stop_ = other.mjd_stop_;
                this->fov_ = other.fov_;
                this->terms_ = other.terms_;
                this->observation_id_ = other.observation_id_;
                this->center_ = other.center_;
                this->meta_ = other.meta_;
            }
            return *this;
        }

        // Property getters
        inline string source() const { return source_; };
        inline string observatory() const { return observatory_; };
        inline string product_id() const { return product_id_; };
        inline optional<int64_t> observation_id() const { return observation_id_; };
        inline double mjd_start() const { return mjd_start_; };
        inline double mjd_stop() const { return mjd_stop_; };
        inline string fov() const { return string(fov_); };
        inline optional<string> center() const { return center_; };
        inline vector<string> terms() const { return terms_; };
        inline optional<string> meta() const { return meta_; };

        // Property setters
        inline void source(const string new_source) { source_ = new_source; };
        inline void observatory(const string name) { observatory_ = name; };
        inline void product_id(const string new_product_id) { product_id_ = new_product_id; };
        void observation_id(optional<int64_t> new_observation_id) { observation_id_ = new_observation_id; };
        inline void mjd_start(double new_mjd_start) { mjd_start_ = new_mjd_start; };
        inline void mjd_stop(double new_mjd_stop) { mjd_stop_ = new_mjd_stop; };
        inline void fov(string new_fov) { fov_ = new_fov; };
        inline void center(optional<string> new_center) { center_ = new_center; };
        void terms(const vector<string> new_terms) { terms_ = new_terms; };
        void terms(const string new_terms) { terms_ = util::split(new_terms, ' '); };
        void meta(const optional<string> new_meta) { meta_ = new_meta; };

        // Calculated properties.

        // Exposure time (s).
        inline double exposure() const { return (mjd_stop_ - mjd_start_) * 86400; };

        // Observation mid-time.
        inline double mjd_mid() const { return (mjd_start_ + mjd_stop_) / 2; };

        // check if observation is valid
        bool is_valid() const;

        // output
        //
        // Show FOV or meta in output?
        struct Format
        {
            bool show_fov = false;
            bool show_meta = false;
        } format;

        friend std::ostream &operator<<(std::ostream &os, const Observation &observation);

        // test if observation has the same FOV as another
        bool is_same_fov(const Observation &other) const;

        // test if observation is equal to another by comparing
        // - source
        // - observatory
        // - product_id
        // - mjd_start
        // - mjd_stop
        // - fov (via is_same_fov)
        // - meta
        bool operator==(const Observation &other) const;
        bool operator!=(const Observation &other) const { return !((*this) == other); };

        // Create a polygon from this observation's field of view, with optional
        // validation checks.
        void as_polygon(S2Polygon &polygon, const bool verify = false) const;

        // Convert to boost JSON object
        json::object as_json();

    private:
        string source_, observatory_, product_id_;
        optional<int64_t> observation_id_;
        double mjd_start_ = 0, mjd_stop_ = 0;
        string fov_;
        optional<string> center_, meta_;
        vector<string> terms_;
    };

    class Observations
    {
    public:
        vector<Observation> data;
        Observation::Format format;

        // Default constructor is an empty vector.
        Observations() {};

        // Initialize with a single Observation
        Observations(const Observation &observation)
        {
            append(observation);
        }

        // Initialize with a vector of Observation
        Observations(const vector<Observation> &observations)
        {
            append(observations);
        }

        // Copy constructor.
        Observations(const Observations &observations)
        {
            append(observations.data);
        };

        // Access element by index.
        Observation &operator[](int i) { return data[i]; };

        // Append a single observation.
        inline void append(const Observation &observation)
        {
            data.push_back(observation);
        };

        // Append a vector of observations.
        inline void append(const vector<Observation> &observations)
        {
            data.reserve(data.size() + observations.size());
            data.insert(data.end(), observations.begin(), observations.end());
        };

        // Append another Observations object.
        inline void append(const Observations &observations)
        {
            append(observations.data);
        }

        // Pointer to beginning of vector.
        auto begin() { return data.begin(); };
        auto begin() const { return data.begin(); };

        // Pointer to end of vector.
        auto end() { return data.end(); };
        auto end() const { return data.end(); };

        // Number of items.
        size_t size() const { return data.size(); }

        // Get all observation_ids.
        vector<optional<int64_t>> observation_ids() const
        {
            vector<optional<int64_t>> ids(size());
            std::transform(begin(), end(), ids.begin(),
                           [](auto const &observation)
                           { return observation.observation_id(); });
            return ids;
        };

        // Remove any duplicate observation IDs, in place.
        void remove_duplicate_observation_ids();
    };

    // Print a table of observations.
    std::ostream &operator<<(std::ostream &os, const Observations &v);
}

// custom specialization of std::hash for unordered_set<Observation>
template <>
struct std::hash<sbsearch::Observation>
{
    std::size_t operator()(sbsearch::Observation const &observation) const noexcept
    {
        return std::hash<int64_t>{}(observation.observation_id().value_or(-1));
    }
};

#endif // OBSERVATION_H_