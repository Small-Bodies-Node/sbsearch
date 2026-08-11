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
        Observation(string_view source,
                    string_view observatory,
                    string_view product_id,
                    const double mjd_start,
                    const double mjd_stop,
                    string_view fov,
                    const vector<string> &terms = {},
                    const optional<int64_t> &observation_id = {},
                    const optional<string> &center = {},
                    const optional<string> &meta = {},
                    const double mjd_added = {});

        // Copy constructor
        Observation(const Observation &other) = default;

        // Move constructor
        Observation(Observation &&other) = default;

        // Copy assignment
        Observation &operator=(const Observation &other) = default;

        // Move assignment
        Observation &operator=(Observation &&other) = default;

        // Property getters
        inline string_view source() const { return source_; };
        inline string_view observatory() const { return observatory_; };
        inline string_view product_id() const { return product_id_; };
        inline optional<int64_t> observation_id() const { return observation_id_; };
        inline double mjd_start() const { return mjd_start_; };
        inline double mjd_stop() const { return mjd_stop_; };
        inline string_view fov() const { return string(fov_); };
        inline optional<string> center() const { return center_; };
        inline vector<string> terms() const { return terms_; };
        inline optional<string> meta() const { return meta_; };
        inline double mjd_added() const { return mjd_added_; };

        // Property setters
        inline void source(const string_view new_source) { source_ = new_source; };
        inline void observatory(const string_view name) { observatory_ = name; };
        inline void product_id(const string_view new_product_id) { product_id_ = new_product_id; };
        void observation_id(optional<int64_t> new_observation_id) { observation_id_ = new_observation_id; };
        inline void mjd_start(double new_mjd_start) { mjd_start_ = new_mjd_start; };
        inline void mjd_stop(double new_mjd_stop) { mjd_stop_ = new_mjd_stop; };
        inline void fov(string_view new_fov) { fov_ = new_fov; };
        inline void center(optional<string> new_center) { center_ = new_center; };
        void terms(const vector<string> new_terms) { terms_ = new_terms; };
        void terms(const string_view new_terms) { terms_ = util::split(new_terms, ' '); };
        void meta(const optional<string> new_meta) { meta_ = new_meta; };
        void mjd_added(const double new_mjd_added) { mjd_added_ = new_mjd_added; };

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
        double mjd_start_ = 0, mjd_stop_ = 0, mjd_added_ = 0;
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
        Observations(Observation &&observation)
        {
            append(std::move(observation));
        }

        // Initialize with a vector of Observation
        Observations(const vector<Observation> &observations)
        {
            append(observations);
        }
        Observations(vector<Observation> &&observations)
        {
            append(std::move(observations));
        }

        // Copy constructor.
        Observations(const Observations &observations) = default;

        // Move constructor.
        Observations(Observations &&observations) = default;

        // Copy assignment.
        Observations &operator=(const Observations &observations) = default;

        // Move assignment.
        Observations &operator=(Observations &&observations) = default;

        // Access element by index.
        Observation &operator[](int i) { return data[i]; };

        // Append a single observation.
        inline void append(const Observation &observation)
        {
            data.push_back(observation);
        };
        inline void append(Observation &&observation)
        {
            data.push_back(std::move(observation));
        };

        // Append a vector of observations.
        inline void append(const vector<Observation> &observations)
        {
            data.reserve(data.size() + observations.size());
            std::copy(observations.begin(), observations.end(), std::back_inserter(data));
        };
        inline void append(vector<Observation> &&observations)
        {
            data.reserve(data.size() + observations.size());
            std::move(observations.begin(), observations.end(), std::back_inserter(data));
        };

        // Append another Observations's object.
        inline void append(const Observations &observations)
        {
            append(observations.data);
        }
        inline void append(Observations &&observations)
        {
            append(std::move(observations.data));
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