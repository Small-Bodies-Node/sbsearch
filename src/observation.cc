#include "config.h"

#include <algorithm>
#include <cinttypes>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <boost/json.hpp>
#include <s2/s2error.h>
#include <s2/s2polygon.h>
#include <s2/s2latlng.h>

#include "observation.h"
#include "exceptions.h"
#include "table.h"
#include "util.h"

using sbsearch::table::Table;
using std::string;
using std::vector;
namespace json = boost::json;

namespace sbsearch
{
    Observation::Observation(string source,
                             string observatory,
                             string product_id,
                             double mjd_start,
                             double mjd_stop,
                             string fov,
                             vector<string> terms,
                             optional<int64_t> observation_id,
                             optional<string> center)
    {
        source_ = source;
        observatory_ = observatory;
        product_id_ = product_id;
        observation_id_ = observation_id;
        mjd_start_ = mjd_start;
        mjd_stop_ = mjd_stop;
        fov_ = string(fov);
        terms_ = terms;
        center_ = center;
        is_valid();
    }

    void Observation::observation_id(optional<int64_t> new_observation_id)
    {
        // To help prevent database corruption, observation IDs can be "erased"
        // but cannot be simply replaced with a new ID.
        if (observation_id_.has_value() & (new_observation_id != std::nullopt))
            throw ObservationError("ID already defined.");
        else
            observation_id_ = new_observation_id;
    };

    bool Observation::is_valid() const
    {
        // checks `fov`: must be parsable into at least 3 vertices
        // ensures that stop >= start
        vector<S2Point> vertices = sbsearch::make_vertices(string(fov_));
        if (vertices.size() < 3)
            throw std::runtime_error("FOV must be parsable into at least three vertices.");

        if (mjd_stop_ < mjd_start_)
            throw std::runtime_error("Observation stops before it starts.");

        return true;
    }

    std::ostream &operator<<(std::ostream &os, const Observation &observation)
    {
        os
            << observation.observation_id().value_or(-1) << "  "
            << '"' << observation.source() << '"'
            << "  "
            << '"' << observation.observatory() << '"'
            << "  "
            << '"' << observation.product_id() << '"'
            << "  "
            << observation.mjd_start() << "  "
            << observation.mjd_stop() << "  "
            << (observation.mjd_stop() - observation.mjd_start()) * 86400;

        if (observation.format.show_fov)
            os << "  "
               << '"' << observation.fov() << '"';

        return os;
    }

    bool Observation::is_same_fov(const Observation &other) const
    {
        S2Polygon this_polygon, other_polygon;
        other.as_polygon(other_polygon);
        as_polygon(this_polygon);
        return this_polygon.BoundaryEquals(other_polygon);
    }

    bool Observation::operator==(const Observation &other) const
    {
        return (
            (source_ == other.source()) &
            (observatory_ == other.observatory()) &
            (product_id_ == other.product_id()) &
            (observation_id_ == other.observation_id()) &
            (mjd_start_ == other.mjd_start()) &
            (mjd_stop_ == other.mjd_stop()) &
            is_same_fov(other));
    }

    void Observation::terms(string new_terms)
    {
        terms_ = split(new_terms, ' ');
    }

    void Observation::terms(vector<string> new_terms)
    {
        terms_ = new_terms;
    }

    void Observation::as_polygon(S2Polygon &polygon) const
    {
        make_polygon(string(fov_), polygon);
    };

    json::object Observation::as_json()
    {
        json::object obj;
        obj["source"] = source();
        obj["observatory"] = observatory();
        obj["product_id"] = product_id();
        if (observation_id_)
            obj["observation_id"] = observation_id_.value();
        else
            obj["observation_id"] = nullptr;
        obj["mjd_start"] = mjd_start();
        obj["mjd_stop"] = mjd_stop();
        obj["fov"] = fov();
        return obj;
    }

    std::ostream &operator<<(std::ostream &os, const Observations &observations)
    {
        int n = observations.size();

        bool show_fov = false;
        vector<string> sources(n), observatories(n), product_ids(n), fovs(n);
        vector<int64_t> observation_ids(n);
        vector<double> mjd_starts(n), mjd_stops(n), exposures(n);

        if (n > 0)
            show_fov = std::max_element(observations.begin(), observations.end(),
                                        [](const Observation &a, const Observation &b)
                                        { return a.format.show_fov < b.format.show_fov; })
                           ->format.show_fov;

        std::transform(observations.begin(), observations.end(), sources.begin(),
                       [](const Observation &obs)
                       { return obs.source(); });

        std::transform(observations.begin(), observations.end(), observatories.begin(),
                       [](const Observation &obs)
                       { return obs.observatory(); });

        std::transform(observations.begin(), observations.end(), product_ids.begin(),
                       [](const Observation &obs)
                       { return obs.product_id(); });

        std::transform(observations.begin(), observations.end(), fovs.begin(),
                       [](const Observation &obs)
                       { return obs.fov(); });

        std::transform(observations.begin(), observations.end(), observation_ids.begin(),
                       [](const Observation &obs)
                       { return obs.observation_id().value_or(-1); });

        std::transform(observations.begin(), observations.end(), mjd_starts.begin(),
                       [](const Observation &obs)
                       { return obs.mjd_start(); });

        std::transform(observations.begin(), observations.end(), mjd_stops.begin(),
                       [](const Observation &obs)
                       { return obs.mjd_stop(); });

        std::transform(observations.begin(), observations.end(), exposures.begin(),
                       [](const Observation &obs)
                       { return (obs.mjd_stop() - obs.mjd_start()) * 86400; });

        Table table;
        table.add_column("observation_id", "%" PRId64, observation_ids);
        table.add_column("source", "%s", sources);
        table.add_column("product_id", "%s", product_ids);
        table.add_column("observatory", "%s", observatories);
        table.add_column("mjd_start", "%.6lf", mjd_starts);
        table.add_column("mjd_stop", "%.6lf", mjd_stops);
        table.add_column("exposure", "%.3lf", exposures);
        if (show_fov)
            table.add_column("fov", "%s", fovs);

        os << table;
        return os;
    }

    void Observations::remove_duplicate_observation_ids()
    {
        std::unordered_set<int64_t> unique_ids;
        unique_ids.reserve(data.size());

        auto is_duplicate = [&](auto const &observation)
        { return !unique_ids.insert(observation.observation_id().value_or(-1)).second; };

        data.erase(std::remove_if(data.begin(), data.end(), is_duplicate), data.end());
    }
}
