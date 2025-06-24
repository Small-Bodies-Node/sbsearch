#include "config.h"

#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>
#include <boost/json.hpp>
#include <s2/s2cell_id.h>
#include <s2/s2cell.h>
#include <s2/s2latlng.h>
#include <s2/s2polygon.h>
#include <s2/s2point.h>

#include "ephemeris.h"
#include "found.h"
#include "query_info.h"
#include "observation.h"
#include "util/string.h"

using std::set;
using std::string;
using std::vector;

namespace sbsearch
{
    QueryInfo::Polygon::Polygon(const vector<S2Point> &vertices)
    {
        assert(vertices.size() == 4);
        for (int i = 0; i < 4; i++)
        {
            const S2LatLng ll(vertices[i]);
            this->at(i) = {ll.lng().degrees(), ll.lat().degrees()};
        }
    }

    QueryInfo::Polygon::Polygon(const std::unique_ptr<S2Polygon> &polygon)
    {
        assert(polygon->loop(0)->num_vertices() == 4);

        vector<S2Point> vertices;
        for (int i = 0; i < 4; i++)
            vertices.emplace_back(polygon->loop(0)->vertex(i));

        *this = Polygon(vertices);
    }

    boost::json::array QueryInfo::Polygon::as_json()
    {
        boost::json::array vertices;
        for (int i = 0; i < 4; i++)
        {
            boost::json::array coordinates = {(*this)[i][0], (*this)[i][1]};
            vertices.emplace_back(coordinates);
        }
        return vertices;
    }

    QueryInfo::Polygons::Polygons(const vector<string> &fovs)
    {
        for (auto const &fov : fovs)
            emplace_back(Polygon(util::make_vertices(fov)));
    }

    QueryInfo::QueryInfo()
    {
        using boost::json::array;
        using boost::json::object;
        data.emplace_object() = {
            {"observations", {{"polygons", object()}, {"terms", object()}}},
            {"matches", array()},
            {"ephemeris", {{"polygons", array()}, {"segments", array()}, {"terms", object()}}}};
    }

    void QueryInfo::approximate_matches(const Observations &observations)
    {
        for (auto const &observation : observations)
        {
            auto const obsid = observation.observation_id();
            if (!obsid)
                continue;

            auto &polygons = data.at_pointer("/observations/polygons").as_object();
            polygons[std::to_string(obsid.value())].emplace_array() = Polygon(observation.fov()).as_json();

            save_terms(observation.terms(), data.at_pointer("/observations/terms").as_object());
        }
    }

    void QueryInfo::matches(const Founds &founds)
    {
        for (auto const &found : founds)
        {
            auto const obsid = found.observation.observation_id();
            if (!obsid)
                continue;

            data.at("matches").as_array().emplace_back(obsid.value());
        }
    }

    void QueryInfo::ephemeris_segment(const Ephemeris &segment,
                                      const double padding,
                                      const vector<string> &query_terms)
    {
        save_terms(query_terms, data.at_pointer("/ephemeris/terms").as_object());

        auto &polygons = data.at("ephemeris").at("polygons").as_array();
        for (auto const &polygon : segment.as_polygons(padding))
            polygons.emplace_back(Polygon(polygon).as_json());

        // save ephemeris positions and uncertainty ellipse
        boost::json::array eph_data;
        for (auto const &point : segment.data())
        {
            boost::json::object eph;
            eph["ra"] = point.ra.value_or(1e99);
            eph["dec"] = point.dec.value_or(1e99);
            eph["unc a"] = point.unc_a.value_or(1e99);
            eph["unc b"] = point.unc_b.value_or(1e99);
            eph["unc theta"] = point.unc_theta.value_or(1e99);
            eph_data.emplace_back(eph);
        }
        data.at("ephemeris").at("segments").as_array().emplace_back(std::move(eph_data));
    }

    void QueryInfo::save_terms(const vector<string> &terms, boost::json::object &dest)
    {
        for (auto const &term : terms)
        {
            if (dest.contains(term))
                continue;

            const int start = (term[0] == '$') ? 1 : 0;
            S2Cell cell(S2CellId::FromToken(term.substr(start)));

            vector<S2Point> vertices;
            for (int i = 0; i < 4; i++)
                vertices.emplace_back(S2LatLng(cell.GetVertex(i)).Normalized());

            auto &array = dest[term].emplace_array();
            array = Polygon(vertices).as_json();
        }
    }

    std::ostream &operator<<(std::ostream &os, const QueryInfo &info)
    {
        return os << info.data;
    }
}
