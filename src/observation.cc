#include "config.h"

#include <algorithm>
#include <ostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <boost/json.hpp>
#include <s2/s2error.h>
#include <s2/s2polygon.h>
#include <s2/s2latlng.h>

#include "observation.h"
#include "table.h"
#include "util.h"

using sbsearch::table::Table;
using std::string;
using std::vector;
namespace json = boost::json;

namespace sbsearch
{
    Observation::Observation(string source, string observatory, string product_id, double mjd_start, double mjd_stop, string fov, string terms, int64 observation_id)
    {
        source_ = source;
        observatory_ = observatory;
        product_id_ = product_id;
        observation_id_ = observation_id;
        mjd_start_ = mjd_start;
        mjd_stop_ = mjd_stop;
        fov_ = string(fov);
        terms_ = string(terms);
        is_valid();
    }

    void Observation::observation_id(int64 new_observation_id)
    {
        // To help prevent database corruption, observation IDs may not be
        // updated if they are already defined.
        if (observation_id_ != UNDEFINED_OBSID)
            throw std::runtime_error("Observation ID already defined.");
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
            << observation.observation_id() << "  "
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
        terms_ = string(new_terms);
    }

    void Observation::terms(vector<string> new_terms)
    {
        terms(join(new_terms, " "));
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
        obj["observation_id"] = observation_id();
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
        vector<int64> observation_ids(n);
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
                       { return obs.observation_id(); });

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
        table.add_column("observatory", "%s", observatories);
        table.add_column("mjd_start", "%.6lf", mjd_starts);
        table.add_column("mjd_stop", "%.6lf", mjd_stops);
        table.add_column("exposure", "%.3lf", exposures);
        if (show_fov)
            table.add_column("fov", "%s", fovs);

        os << table;
        return os;
    }
}
