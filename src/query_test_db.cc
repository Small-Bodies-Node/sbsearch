#include "test_db.h"
#include "sbsearch.h"

#include <cstdio>
#include <cstdlib>
#include <set>
#include <sstream>
#include <iostream>
#include <list>
#include <vector>
#include <string>
#include <random>
#include <ctime>

#include <sqlite3.h>

#include "s2/s2region_term_indexer.h"
#include "s2/s2polyline.h"
#include "s2/s2metrics.h"
#include "s2/s2latlng.h"
#include "s2/s2builder.h"
#include "s2/s2error.h"
#include "s2/s2testing.h"
#include "s2/s2builderutil_s2polygon_layer.h"
#include "s2/s2text_format.h"

#include "util.h"
#include "observation.h"
#include "ephemeris.h"
#include "sbsdb_sqlite3.h"

#define N_COMETS 1

using sbsearch::Ephemeris;
using sbsearch::mjd_to_time_terms;
using sbsearch::Observation;
using sbsearch::SBSearchDatabaseSqlite3;
using sbsearch::sql_check;
using std::cerr;
using std::cout;
using std::endl;
using std::vector;

struct Found
{
    Observation obs;
    Ephemeris eph;
};

Ephemeris get_ephemeris(const double mjd0, const double mjd1, const double step,
                        const double ra0, const double dec0,
                        const double ra_rate, const double dec_rate)
{
    double ra = ra0, dec = dec0;
    double dec_dir = 1;
    vector<S2Point> all_vertices;
    vector<double> all_times;
    const S1Angle tol = S1Angle::Radians(EPHEMERIS_TOLERANCE);

    for (double mjd = mjd0; mjd <= mjd1; mjd += step)
    {
        all_vertices.push_back(S2LatLng::FromDegrees(dec, ra).ToPoint());
        all_times.push_back(mjd);
        ra += ra_rate * std::cos(dec * PI / 180) * step;
        dec += dec_dir * dec_rate * step;
        if (dec > 90)
        {
            dec = 90 - (dec - 90);
            dec_dir *= -1;
            ra += 180;
        }
        else if (dec < -90)
        {
            dec = -90 - (dec + 90);
            dec_dir *= -1;
            ra += 180;
        }

        if (ra > 180)
        {
            ra -= 360;
        }
        else if (ra < -180)
        {
            ra += 360;
        }
    }
    S2Polyline full_resolution(all_vertices);

    // simplify
    vector<int> indices;
    vector<S2Point> vertices;
    vector<double> times;
    full_resolution.SubsampleVertices(tol, &indices);
    for (const auto &i : indices)
    {
        vertices.push_back(all_vertices[i]);
        times.push_back(all_times[i]);
    }

    return Ephemeris(vertices, times);
}

/*
static int collect_found_rowids(void *found_ptr, int count, char **data, char **columns)
{
    static int total = 0;
    total += count;

    std::set<int64> *found = static_cast<std::set<int64> *>(found_ptr);

    for (int i = 0; i < count; i++)
        found->insert(std::strtoll(data[i], NULL, 10));

    return 0;
}

std::set<int64> fuzzy_search(sqlite3 *db, S2RegionTermIndexer &indexer, Ephemeris eph)
{
    // Collect search terms
    std::set<string> terms;
    for (auto segment : eph.segments())
    {
        for (auto term : segment.query_terms(indexer))
            terms.insert(term);
    }

    char *error_message = 0;
    char statement[MAXIMUM_QUERY_CLAUSE_LENGTH + 100];
    char *end = statement;
    int maximum_term_string_length;
    int count = 0;
    string term_string;
    std::set<int64> approximate_matches;

    // Query database with terms, but not too many at once
    end = stpcpy(statement, "SELECT rowid FROM obs_geometry_time WHERE terms MATCH '");
    for (auto term : terms)
    {
        if (++count % 100 == 0)
            cout << "\r  - " << count << " of " << terms.size() << " query terms." << std::flush;

        if (end != (statement + 55))
        {
            end = stpcpy(end, " OR ");
        }
        end = stpcpy(end, term.c_str());

        if ((end - statement) > MAXIMUM_QUERY_CLAUSE_LENGTH)
        {
            strcpy(end, "';");
            sql_check(sqlite3_exec(db, statement, collect_found_rowids, &approximate_matches, &error_message), error_message);
            end = statement + 55;
        }
    }
    cout << "\r  - " << count << " of " << terms.size() << " query terms." << endl;
    return approximate_matches;
}

vector<Found> find_intersecting_observations(sqlite3 *db, std::set<int64> obsids, Ephemeris eph)
{
    vector<Found> found;
    S2Polyline polyline = eph.as_polyline();
    for (auto obsid : obsids)
    {
        Observation obs(db, obsid);
        unique_ptr<S2Polygon> polygon = obs.as_polygon();
        if (polygon->Intersects(polyline))
            found.push_back(Found{obs, eph});
    }
    return found;
}
*/

vector<Found> query_ephemeris(SBSearchDatabaseSqlite3 &db, Ephemeris eph)
{
    assert(eph.num_segments() > 0);
    cout << "\n     - Querying " << eph.num_segments() << " ephemeris segments." << endl;

    vector<Observation> approximate_matches = db.fuzzy_search(eph);
    cout << "     - Fuzzy search found " << approximate_matches.size() << " observations." << endl;
    // vector<Found> found = find_intersecting_observations(db, approximate_matches, eph);
    // cout << "\nIntersection test found " << found.size() << " observations." << endl;
    return vector<Found>();
}

std::pair<double, double> survey_date_range(SBSearchDatabaseSqlite3 &db)
{
    double mjd_start, mjd_stop;
    mjd_start = db.get_one_value<double>("SELECT MIN(mjd_start) FROM observations;");
    mjd_stop = db.get_one_value<double>("SELECT MAX(mjd_stop) FROM observations;");

    return std::pair<double, double>(mjd_start, mjd_stop);
}

void query_test_db()
{
    // get date range for query
    SBSearchDatabaseSqlite3 db("sbsearch_test.db");
    std::pair<double, double> date_range = survey_date_range(db);
    cout << "\nGenerating " << (date_range.second - date_range.first) / 365.25 << " year long ephemerides:\n";

    // build ephemeris
    double step = 1.0; // days
    double ra0 = FOV_WIDTH * 1.5;
    double dec0 = 0;
    double ra_rate, dec_rate;
    std::srand((unsigned)std::time(0));
    for (int i = 0; i < N_COMETS; ++i)
    {                                                                       // -100 to 100 arcsec/hr
        dec_rate = ((float)std::rand() / RAND_MAX - 0.5) * 200 / 3600 * 24; // deg/day
        ra_rate = ((float)std::rand() / RAND_MAX - 0.5) * 200 / 3600 * 24;
        printf("%3d: μ = %lf deg/day", i + 1, std::hypot(ra_rate, dec_rate));

        Ephemeris eph = get_ephemeris(date_range.first, date_range.second, step, ra0, dec0, ra_rate, dec_rate);
        vector<Found> found = query_ephemeris(db, eph);
    }
    cout << "\n\n";
}

int main(int argc, char **argv)
{
    /*
    int rc = 0;
    sqlite3 *db;
    vector<string> terms;

    S2RegionTermIndexer::Options options;
    options.set_max_level(S2::kAvgEdge.GetClosestLevel(3e-4));
    S2RegionTermIndexer indexer(options);

    rc = sqlite3_open("sbsearch_test.db", &db);
    if (rc != SQLITE_OK)
    {
        cerr << "Error opening database" << sqlite3_errmsg(db) << endl;
        return (-1);
    }
    else
        cout << "Opened Database Successfully!" << endl;

    std::pair<double, double> date_range = survey_date_range(db);
    cout << "Generating " << (date_range.second - date_range.first) / 365.25 << " year ephemeris." << endl;

    const double step = 1.0; // days
    double ra0 = FOV_WIDTH * 1.5;
    double dec0 = 0;
    double ra_rate, dec_rate;
    std::srand((unsigned)std::time(0));
    for (int i = 0; i < N_COMETS; ++i)
    {
        // S2LatLng start(S2Testing::RandomPoint());
        // -100 to 100 arcsec/hr
        dec_rate = ((float)std::rand() / RAND_MAX - 0.5) * 200 / 3600 * 24; // deg/day
        ra_rate = ((float)std::rand() / RAND_MAX - 0.5) * 200 / 3600 * 24;

        cout << "Generating ephemeris\n";
        cout << "  dRA = " << ra_rate << "deg/day\n  dDec = " << dec_rate << "deg/day\n";
        Ephemeris eph = get_ephemeris(date_range.first, date_range.second, step, ra0, dec0, ra_rate, dec_rate);
        vector<Found> found = query_ephemeris(db, indexer, eph);

        std::for_each_n(found.begin(), std::min<int>(found.size(), 30), [](auto &f)
                        { cout << f.obs.obsid() << " " << f.obs.mjd_start() << " " << f.obs.mjd_stop() << " " << f.obs.fov() << endl; });
        if (found.size() > 30)
            cout << "...\n";
    }

    sqlite3_close(db);
    */
    query_test_db();
    return 0;
}