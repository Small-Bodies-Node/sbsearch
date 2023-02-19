#include "config.h"

#include <algorithm>
#include <ostream>
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
    Observation::Observation(string source, string product_id, double mjd_start, double mjd_stop, string fov, string terms, int64 observation_id)
    {
        source_ = source;
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
        vector<S2Point> vertices = sbsearch::makeVertices(string(fov_));
        if (vertices.size() < 3)
            throw std::runtime_error("FOV must be parsable into at least three vertices.");

        if (mjd_stop_ < mjd_start_)
            throw std::runtime_error("Observation stops before it starts.");

        return true;
    }

    Observation::Format Observation::format_widths() const
    {
        double exposure_time = (mjd_stop() - mjd_start()) * 86400;
        Observation::Format format = {
            size_t(std::floor(std::log10(observation_id()))) + 1,
            source().length(),
            product_id().length(),
            size_t(std::floor(std::log10(exposure_time))) + 3,
            fov().length(),
            format.show_fov,
            format.quote_strings};
        return format;
    }

    std::ostream &operator<<(std::ostream &os, const Observation &observation)
    {
        os << std::right
           << std::fixed
           << std::setw(observation.format.observation_id_width)
           << observation.observation_id() << "  "
           << (observation.format.quote_strings ? "\"" : "")
           << std::setw(observation.format.source_width)
           << observation.source()
           << (observation.format.quote_strings ? "\"" : "")
           << "  "
           << (observation.format.quote_strings ? "\"" : "")
           << std::setw(observation.format.product_id_width)
           << observation.product_id()
           << (observation.format.quote_strings ? "\"" : "")
           << "  "
           << std::setw(11)
           << std::setprecision(5)
           << observation.mjd_start() << "  "
           << std::setw(11)
           << std::setprecision(5)
           << observation.mjd_stop() << "  "
           << std::setw(observation.format.exposure_time_width)
           << std::setprecision(1)
           << (observation.mjd_stop() - observation.mjd_start()) * 86400;

        if (observation.format.show_fov)
            os << "  "
               << (observation.format.quote_strings ? "\"" : "")
               << std::setw(observation.format.fov_width)
               << observation.fov()
               << (observation.format.quote_strings ? "\"" : "");

        os << std::defaultfloat;

        return os;
    }

    bool Observation::is_same_fov(const Observation &other) const
    {
        auto other_polygon = other.as_polygon();
        return as_polygon().BoundaryEquals(other_polygon);
    }

    bool Observation::is_equal(const Observation &other) const
    {
        return ((source_ == other.source()) &
                (product_id_ == other.product_id()) &
                is_same_fov(other) & (mjd_start_ == other.mjd_start()) &
                (mjd_stop_ == other.mjd_stop()) &
                (observation_id_ == other.observation_id()));
    }

    bool Observation::operator==(const Observation &other) const
    {
        return is_equal(other);
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

    std::ostream &operator<<(std::ostream &os, const vector<Observation> &observations)
    {
        // scan vector to determine column widths
        Observation::Format format;
        for (const Observation &observation : observations)
        {
            Observation::Format _format = observation.format_widths();
            format.observation_id_width = std::max(format.observation_id_width, _format.observation_id_width);
            format.source_width = std::max(format.source_width, _format.source_width);
            format.product_id_width = std::max(format.product_id_width, _format.product_id_width);
            format.fov_width = std::max(format.fov_width, _format.fov_width);
        }

        for (Observation observation : observations)
        {
            observation.format = format;
            os << observation << "\n";
        }
        return os;
    }
}