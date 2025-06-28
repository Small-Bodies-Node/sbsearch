#ifndef QUERYINFO_H_
#define QUERYINFO_H_

#include <array>
#include <mutex>
#include <set>
#include <string>
#include <tuple>
#include <vector>
#include <boost/json.hpp>
#include <s2/s2polygon.h>
#include <s2/s2point.h>

#include "ephemeris.h"
#include "found.h"
#include "observation.h"
#include "util/string.h"

using std::array;
using std::set;
using std::string;
using std::vector;

namespace sbsearch
{
    // Information that may be used to debug or visualize an SBSearch query.
    //
    // Let data be the JSON object, then:
    //
    // * data["observations"] : Observations returned from approximate search.
    //   * ["polygons"] : The observation fields of view, keyed by observation
    //     ID.
    //   * ["terms"] : The observation cells keyed by index term.
    //
    // * data["matches"] : Observation IDs that matched the query.
    //
    // * data["ephemeris"] : Target ephemeris data.
    //   * ["polygons"] : Array of ephemeris query areas.
    //   * ["segments"] : Array of ephemeris segments used in the search: ra,
    //     dec, unc a, unc b, unc theta (deg, deg, arcsec, arcsec, deg).
    //   * ["terms"] : Ephemeris cells keyed by query term.
    //
    // Polygons are (ra, dec) pairs in degrees, 0 to 360.
    class QueryInfo
    {
    public:
        typedef array<double, 2> Coordinates;

        struct Polygon : public array<Coordinates, 4>
        {
            Polygon(const vector<S2Point> &vertices);
            Polygon(const std::unique_ptr<S2Polygon> &polygon);
            Polygon(std::string_view fov) : Polygon(util::make_vertices(fov)) {};
            boost::json::array as_json();
        };

        struct Polygons : public vector<Polygon>
        {
            Polygons(const vector<string> &fovs);
        };

        QueryInfo();

        // Move assignment.  The access mutex is not moved.
        QueryInfo &operator=(QueryInfo &&other)
        {
            if (this != &other)
                this->data = std::move(other.data);

            return *this;
        }

        boost::json::value data;

        // Save observation polygons and terms to data["observations"].
        void approximate_matches(const Observations &observations);

        // Save found observation IDs to data["matches"].
        void matches(const Founds &founds);

        // Save ephemeris segment data to data["ephemeris"].
        void ephemeris_segment(const Ephemeris &ephemeris,
                               const double padding,
                               const vector<string> &query_terms);

    private:
        std::mutex access;
        void save_terms(const vector<string> &terms, boost::json::object &dest);
    };

    std::ostream &operator<<(std::ostream &os, const sbsearch::QueryInfo &info);
}

#endif