#include "config.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <s2/s2latlng.h>
#include <s2/s2loop.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>

#include "logging.h"
#include "string.h"

using std::string;
using std::vector;

namespace sbsearch::util
{
    vector<string> split(std::string_view str, const char delimiter)
    {
        int start = 0, end;
        vector<string> parts;
        while ((end = str.find(delimiter, start)) != string::npos)
        {
            parts.emplace_back(str.substr(start, end - start));
            start = end + 1;
        }

        if (start < str.length())
            parts.emplace_back(str.substr(start)); // remainder of string

        return parts;
    }

    string strip(std::string_view str)
    {
        string::size_type start = str.find_first_not_of(" ");
        string::size_type stop = str.find_last_not_of(" ");
        if (stop == string::npos)
            return "";
        return string(str.substr(start, stop - start + 1));
    }

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
        for (const string &coord : split(fov, ','))
        {
            vector<string> values = split(coord, ':');
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