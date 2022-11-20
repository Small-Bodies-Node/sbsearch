#include "observation.h"
#include "util.h"

#include <string>
#include <numeric>
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
    Observation::Observation(int64 obsid, double mjd_start, double mjd_stop, const char *fov)
    {
        obsid_ = obsid;
        mjd_start_ = mjd_start;
        mjd_stop_ = mjd_stop;
        strncpy(fov_, fov, 511);
        is_valid();
    }

    Observation::Observation(int64 obsid, double mjd_start, double mjd_stop, S2LatLngRect fov)
    {
        obsid_ = obsid;
        mjd_start_ = mjd_start;
        mjd_stop_ = mjd_stop;

        // field of view as set of comma-separated RA:Dec pairs in degrees
        sprintf(fov_, "%f:%f, %f:%f, %f:%f, %f:%f",
                fov.lat_lo().degrees(), fov.lng_lo().degrees(),
                fov.lat_lo().degrees(), fov.lng_hi().degrees(),
                fov.lat_hi().degrees(), fov.lng_hi().degrees(),
                fov.lat_hi().degrees(), fov.lng_lo().degrees());
    }

    Observation::Observation(sqlite3 *db, int64 obsid)
    {
        char *error_message = 0;
        sqlite3_stmt *statement;

        sqlite3_prepare_v2(db, "SELECT obsid, mjdstart, mjdstop, fov FROM obs WHERE obsid = ?;", -1, &statement, NULL);
        int rc = sqlite3_bind_int64(statement, 1, obsid);

        if (rc != SQLITE_OK)
        {
            std::cerr << sqlite3_errmsg(db) << std::endl;
            throw "Error preparing SQL statement";
        }

        rc = sqlite3_step(statement);
        if (rc == SQLITE_DONE)
        {
            std::cerr << "obsid " << obsid;
            throw "No matching rows for obsid.";
        }

        if (rc != SQLITE_ROW)
        {
            std::cerr << sqlite3_errmsg(db) << std::endl;
            throw "Error retrieving observation";
        }

        obsid_ = sqlite3_column_int64(statement, 0);
        mjd_start_ = sqlite3_column_double(statement, 1);
        mjd_stop_ = sqlite3_column_double(statement, 2);
        strncpy(fov_, (char *)sqlite3_column_text(statement, 3), 511);

        sqlite3_finalize(statement);
    };

    bool Observation::is_valid()
    {
        // presently checks `fov`: must be parsable into at least 3 vertices
        vector<S2Point> vertices = sbsearch::makeVertices(string(fov_));
        if (vertices.size() < 3)
        {
            throw "FOV must be parsable into at least three vertices.";
        }
        return true;
    }

    vector<string> Observation::index_terms(S2RegionTermIndexer &indexer)
    {
        vector<string> terms;
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

        S2LatLngRect fov = S2LatLngRect::FromPointPair(
            S2LatLng::FromDegrees(corners[0], corners[1]),
            S2LatLng::FromDegrees(corners[4], corners[5]));
        vector<string> spatial_terms = indexer.GetIndexTerms(fov, "");

        // Get terms for the time
        vector<string> time_terms = mjd_to_time_terms(mjd_start_, mjd_stop_);

        // Join query terms, each segment gets a time suffix
        for (auto time_term : time_terms)
            for (auto spatial_term : spatial_terms)
                terms.push_back(spatial_term + "-" + time_term);

        return terms;
    }

    void Observation::add_sql(S2RegionTermIndexer &indexer, char *statement)
    {
        vector<string> terms = index_terms(indexer);

        auto join = [](string a, string b)
        {
            return std::move(a) + " OR " + std::move(b);
        };
        std::string termString = std::accumulate(std::next(terms.begin()), terms.end(), terms[0], join);
        sprintf(statement, "INSERT INTO obs VALUES (%lld, %f, %f, '%s'); INSERT INTO obs_geometry_time(rowid, terms) VALUES (%lld, \"%s\");",
                obsid_, mjd_start_, mjd_stop_, fov_, obsid_, termString.c_str());
    }

    void Observation::add_to_database(sqlite3 *db, S2RegionTermIndexer &indexer)
    {
        char statement[1024];
        add_sql(indexer, statement);
        sql_execute(db, statement);
    }

    unique_ptr<S2Polygon> Observation::as_polygon()
    {
        // std::cerr << fov_ << std::endl;
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
}