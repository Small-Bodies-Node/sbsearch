#include "config.h"

#include <stdexcept>
#include <string>

#include <s2/s2error.h>
#include <s2/s2polygon.h>
#include <s2/s2latlng.h>

#include "observation.h"
#include "util.h"

using sbsearch::format_vertices;
using sbsearch::makePolygon;
using std::string;
using std::vector;

namespace sbsearch
{
    Observation::Observation(double mjd_start, double mjd_stop, string fov, string terms, int64 observation_id)
    {
        observation_id_ = observation_id;
        mjd_start_ = mjd_start;
        mjd_stop_ = mjd_stop;
        fov_ = string(fov);
        terms_ = string(terms);
        is_valid();
    }

    Observation::Observation(double mjd_start, double mjd_stop, vector<S2LatLng> vertices, string terms, int64 observation_id)
    {
        observation_id_ = observation_id;
        mjd_start_ = mjd_start;
        mjd_stop_ = mjd_stop;
        fov_ = format_vertices(vertices);
        terms_ = string(terms);
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
        vector<S2Point> vertices = sbsearch::makeVertices(string(fov_));
        if (vertices.size() < 3)
            throw std::runtime_error("FOV must be parsable into at least three vertices.");

        if (mjd_stop_ < mjd_start_)
            throw std::runtime_error("Observation stops before it starts.");

        return true;
    }

    bool Observation::is_same_fov(Observation &other) const
    {
        auto other_polygon = other.as_polygon();
        return as_polygon().BoundaryEquals(other_polygon);
    }

    bool Observation::is_equal(Observation &other) const
    {
        return (is_same_fov(other) & (mjd_start_ == other.mjd_start()) & (mjd_stop_ == other.mjd_stop()) & (observation_id_ == other.observation_id()));
    }

    void Observation::terms(string new_terms)
    {
        terms_ = string(new_terms);
    }

    void Observation::terms(vector<string> new_terms)
    {
        terms(join(new_terms, " "));
    }

    S2Polygon Observation::as_polygon() const
    {
        S2Polygon polygon;
        makePolygon(string(fov_), polygon);
        return polygon;
    };
}