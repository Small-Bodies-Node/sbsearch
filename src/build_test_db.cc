#include "build_test_db.h"
#include "sbsearch.h"

#include <cstdio>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <numeric>

#include <sqlite3.h>

#include "s2/s2region_term_indexer.h"
#include "s2/s2latlng.h"
#include "s2/s2latlng_rect.h"
#include "s2/s2metrics.h"

#include "observation.h"
#include "util.h"

using sbsearch::mjd_to_time_terms;
using sbsearch::Observation;
using sbsearch::sql_execute;
using std::vector;

int setup_tables(sqlite3 *db)
{

    sql_execute(db, "DROP TABLE IF EXISTS obs;");
    sql_execute(db, "DROP TABLE IF EXISTS obs_geometry_time;");
    sql_execute(db, "CREATE TABLE obs ("
                    "  obsid INTEGER PRIMARY KEY,"
                    "  mjdstart FLOAT NOT NULL,"
                    "  mjdstop FLOAT NOT NULL,"
                    "  fov TEXT NOT NULL"
                    ");");
    sql_execute(db, "CREATE VIRTUAL TABLE obs_geometry_time USING fts4(terms TEXT)");
    sql_execute(db, "DROP INDEX IF EXISTS idx_obs_mjdstart; DROP INDEX IF EXISTS idx_obs_mjdstop;");

    std::cout << "Tables are set." << std::endl;
    return (0);
}

void new_fov(S2LatLngRect &fov)
{
    const static S2LatLng size = S2LatLng::FromDegrees(FOV_WIDTH, FOV_WIDTH);
    const static S1Angle step = S1Angle::Degrees(FOV_WIDTH * 2);
    const static S1Angle two_pi = S1Angle::Radians(2 * PI);
    const static S1Angle pi = S1Angle::Radians(PI);
    static S1Angle ra = S1Angle::Zero();
    static S1Angle dec = S1Angle::Zero();

    ra += step;
    if (ra.radians() > PI)
    {
        dec += step;
        ra -= two_pi;
    }

    if (dec.radians() > PI_2)
    {
        dec -= pi;
    }

    fov = S2LatLngRect::FromCenterSize(S2LatLng(dec, ra).Normalized(), size);
}

int main(int argc, char **argv)
{
    char statement[1024];
    int rc = 0;
    int64 obsid = 0;
    const double mjd0 = 59103.0;
    const char *filename = "sbsearch_test.db";
    double mjd;
    double ra = 0, dec = 0;
    std::vector<string> terms;
    sqlite3 *db;
    S2LatLngRect fov = S2LatLngRect();
    S2RegionTermIndexer::Options options;
    options.set_max_level(S2::kAvgEdge.GetClosestLevel(SPATIAL_TERM_RESOLUTION * 0.00029089));
    options.set_max_cells(8);
    S2RegionTermIndexer indexer(options);

    std::cout << "\nIndex setup:"
              << "\n  Spatial min level: " << options.min_level()
              << "\n  Spatial max level: " << options.max_level()
              << "\n  Time resolution: " << 24.0 / TIME_TERMS_PER_DAY << " hr\n"
              << "\nSurvey setup:"
              << "\n  Nights: " << NIGHTS
              << "\n  Exposures per night: " << EXPOSURES_PER_NIGHT
              << "\n  Cadence: " << CADENCE * 86400 << " s"
              << "\n  Exposure time: " << EXPOSURE * 86400 << " s\n\n";

    rc = sqlite3_open(filename, &db);
    if (rc != SQLITE_OK)
    {
        std::cerr << "Error opening database" << sqlite3_errmsg(db) << std::endl;
        return (-1);
    }
    else
        std::cout << "Opened " << filename << std::endl;
    sql_execute(db, "PRAGMA temp_store_directory = './';");
    setup_tables(db);

    for (int night = 0; night < NIGHTS; night++)
    {
        mjd = mjd0 + night;
        sql_execute(db, "BEGIN TRANSACTION;");
        std::cout << "\rnight " << night + 1 << std::flush;
        for (int i = 0; i < EXPOSURES_PER_NIGHT; ++i, ++obsid)
        {
            mjd += CADENCE;
            new_fov(fov);
            Observation obs(obsid, mjd, mjd + EXPOSURE, fov);
            obs.add_to_database(db, indexer);
        }
        sql_execute(db, "END TRANSACTION;");
    }

    std::cout << std::endl;

    std::cout << "Creating indices\t";
    sql_execute(db, "CREATE INDEX idx_obs_mjdstart ON obs(mjdstart);");
    if (rc != 0)
        return (rc);

    sql_execute(db, "CREATE INDEX idx_obs_mjdstop ON obs(mjdstop);");
    if (rc != 0)
        return (rc);
    std::cout << "Done." << std::endl;

    sqlite3_close(db);
    return (0);
}