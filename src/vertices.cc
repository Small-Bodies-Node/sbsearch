#include <algorithm>
#include <charconv>
#include <string>
#include <string_view>
#include <vector>
#include <s2/s2latlng.h>
#include <s2/s2loop.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>

#include "config.h"
#include "logging.h"
#include "util/string.h"

using std::string;
using std::string_view;
using std::vector;

namespace sbsearch
{
    string format_vertices(const vector<S2LatLng> &vertices)
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

    string format_vertices(const vector<S2Point> &vertices)
    {
        vector<S2LatLng> ll_vertices;
        for (auto vertex : vertices)
            ll_vertices.push_back(S2LatLng(vertex));
        return format_vertices(ll_vertices);
    }

    string format_vertices(const S2LatLngRect &fov)
    {
        vector<S2LatLng> vertices;
        for (int i = 0; i < 4; i++)
            vertices.push_back(fov.GetVertex(i));
        return format_vertices(vertices);
    }

    string format_vertices(const S2Polygon &polygon)
    {
        if (polygon.num_loops() == 0)
            return string();

        vector<S2LatLng> vertices;
        const S2Loop *loop = polygon.loop(0);
        for (int i = 0; i < loop->num_vertices(); i++)
            vertices.push_back(S2LatLng(loop->vertex(i)));
        return format_vertices(vertices);
    }

    vector<S2Point> make_vertices(std::string_view fov)
    {
        vector<S2Point> vertices;
        vertices.reserve(std::count(fov.cbegin(), fov.cend(), ',') + 1);
        for (string_view coord : util::split(fov, ','))
        {
            vector<string> values = util::split(coord, ':');
            if (values.size() < 2)
            {
                Logger::error() << "fov error on vertex \"" << coord << "\" from string \"" << fov << '"' << std::endl;
                throw std::runtime_error("Not enough vertices");
            }

            try
            {
                S2LatLng ll = S2LatLng::FromDegrees(std::stod(values[1]), std::stod(values[0]));
                vertices.push_back(ll.Normalized().ToPoint());
            }
            catch (std::invalid_argument const &err)
            {
                Logger::error() << "fov error on " << coord << " from " << fov << std::endl;
                throw std::runtime_error("Could not parse fov into vertices");
            }
        }
        return vertices;
    }
}