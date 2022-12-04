#include "observation.h"
#include "util.h"
#include "sbsearch.h"

#include <string>
#include <numeric>
#include <stdexcept>
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
    Observation::Observation(double mjd_start, double mjd_stop, const char *fov, string terms, int64 observation_id)
    {
        observation_id_ = observation_id;
        mjd_start_ = mjd_start;
        mjd_stop_ = mjd_stop;
        strncpy(fov_, fov, 511);
        terms_ = string(terms);
        is_valid();
    }

    Observation::Observation(double mjd_start, double mjd_stop, double *lat, double *lng, string terms, int64 observation_id)
    {
        observation_id_ = observation_id;
        mjd_start_ = mjd_start;
        mjd_stop_ = mjd_stop;
        terms_ = string(terms);

        // field of view as set of comma-separated RA:Dec pairs in degrees
        sprintf(fov_, "%f:%f, %f:%f, %f:%f, %f:%f",
                lat[0], lng[0],
                lat[1], lng[1],
                lat[2], lng[2],
                lat[3], lng[3]);
    }

    Observation::Observation(double mjd_start, double mjd_stop, S2LatLngRect fov, string terms, int64 observation_id)
    {
        observation_id_ = observation_id;
        mjd_start_ = mjd_start;
        mjd_stop_ = mjd_stop;
        terms_ = string(terms);

        sprintf(fov_, "%f:%f, %f:%f, %f:%f, %f:%f",
                fov.lat_lo().degrees(), fov.lng_lo().degrees(),
                fov.lat_lo().degrees(), fov.lng_hi().degrees(),
                fov.lat_hi().degrees(), fov.lng_hi().degrees(),
                fov.lat_hi().degrees(), fov.lng_lo().degrees());
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

    void Observation::copy_fov(char *fov)
    {
        for (int i = 0; i < 512; i++)
        {
            fov_[i] = fov[i];
            if (fov[i] == '\0')
                break;
        }
    };

    vector<string> Observation::generate_terms(TermStyle style, S2RegionTermIndexer &indexer)
    {
        vector<double> corners;
        char value[32];
        int start = 0;
        for (int pos = 0, i = 0; pos < 512; pos++, i++)
        {
            if ((fov_[pos] == ':') | (fov_[pos] == ',') | (fov_[pos] == '\0'))
            {
                value[i] = '\0';
                corners.push_back(atof(value));
                if (corners.size() == 8)
                    break;
                i = -1;
            }
            else
            {
                value[i] = fov_[pos];
            }
        }

        // TODO: replace with appropriate polynomial
        S2LatLngRect fov = S2LatLngRect::FromPointPair(
            S2LatLng::FromDegrees(corners[0], corners[1]),
            S2LatLng::FromDegrees(corners[4], corners[5]));
        vector<string> spatial_terms = (style == Observation::index)
                                           ? indexer.GetIndexTerms(fov, "")
                                           : indexer.GetQueryTerms(fov, "");

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