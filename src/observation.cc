#include "observation.h"
#include "util.h"
#include "sbsearch.h"

#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>

#include <sqlite3.h>
#include <s2/s2error.h>
#include <s2/s2polygon.h>
#include <s2/s2builder.h>
#include <s2/s2builderutil_s2polygon_layer.h>
#include <s2/s2latlng.h>
#include <s2/s2latlng_rect.h>
#include <s2/s2region_term_indexer.h>

using sbsearch::makePolygon;
using sbsearch::sql_execute;
using std::string;
using std::unique_ptr;
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
        fov_ = Observation::format_vertices(vertices);
        terms_ = string(terms);
    }

    string Observation::format_vertices(vector<S2LatLng> vertices)
    {
        // field of view as set of comma-separated RA:Dec pairs in degrees
        string fov;
        for (auto vertex : vertices)
        {
            fov += std::to_string(vertex.lng().degrees()) + ":" + std::to_string(vertex.lat().degrees());
            if (vertex != *(vertices.end() - 1))
                fov += ", ";
        }
        return fov;
    }

    string Observation::format_vertices(vector<S2Point> vertices)
    {
        vector<S2LatLng> ll_vertices;
        for (auto vertex : vertices)
            ll_vertices.push_back(S2LatLng(vertex));
        return format_vertices(ll_vertices);
    }

    string Observation::format_vertices(S2LatLngRect fov)
    {
        vector<S2LatLng> vertices;
        for (int i = 0; i < 4; i++)
            vertices.push_back(fov.GetVertex(i));
        return format_vertices(vertices);
    }

    string Observation::format_vertices(int num_vertices, double *ra, double *dec)
    {
        vector<S2LatLng> vertices;
        for (int i = 0; i < num_vertices; i++)
            vertices.push_back(S2LatLng::FromDegrees(dec[i], ra[i]));
        return format_vertices(vertices);
    }

    bool Observation::is_valid()
    {
        // presently checks `fov`: must be parsable into at least 3 vertices
        vector<S2Point> vertices = sbsearch::makeVertices(string(fov_));
        if (vertices.size() < 3)
        {
            throw std::runtime_error("FOV must be parsable into at least three vertices.");
        }
        return true;
    }

    bool Observation::is_same_fov(Observation &other)
    {
        auto other_polygon = other.as_polygon();
        return as_polygon()->BoundaryEquals(*other_polygon.get());
    }

    bool Observation::is_equal(Observation &other)
    {
        return (is_same_fov(other) & (mjd_start_ == other.mjd_start()) & (mjd_stop_ == other.mjd_stop()) & (observation_id_ == other.observation_id()));
    }

    vector<string> Observation::index_terms(S2RegionTermIndexer &indexer, bool update)
    {
        vector<string> terms = generate_terms(Observation::index, indexer);

        if (update)
            terms_ = join(terms, " ");

        return terms;
    }

    vector<string> Observation::query_terms(S2RegionTermIndexer &indexer)
    {
        return generate_terms(Observation::query, indexer);
    }

    unique_ptr<S2Polygon> Observation::as_polygon()
    {
        return makePolygon(string(fov_));
    };

    vector<string> Observation::generate_terms(TermStyle style, S2RegionTermIndexer &indexer)
    {
        auto polygon = as_polygon();
        vector<string> spatial_terms = (style == Observation::index)
                                           ? indexer.GetIndexTerms(*polygon, "")
                                           : indexer.GetQueryTerms(*polygon, "");

        // Get terms for the time
        vector<string> time_terms = mjd_to_time_terms(mjd_start_, mjd_stop_);

        // Join query terms, each segment gets a time suffix, save to terms string
        vector<string> terms;
        for (auto time_term : time_terms)
            for (auto spatial_term : spatial_terms)
                terms.push_back(spatial_term + "-" + time_term);

        return terms;
    }
}