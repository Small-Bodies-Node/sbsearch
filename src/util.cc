#include "util.h"
#include "sbsearch.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <ostream>
#include <memory>
#include <string>
#include <vector>

#include <sqlite3.h>

#include <s2/s1angle.h>
#include <s2/s2builder.h>
#include <s2/s2builderutil_s2polygon_layer.h>
#include <s2/s2latlng.h>
#include <s2/s2loop.h>
#include <s2/s2point.h>

using std::atan2;
using std::ceil;
using std::cos;
using std::floor;
using std::sin;
using std::string;
using std::tan;
using std::vector;

namespace sbsearch
{
    vector<string> mjd_to_time_terms(const double start, const double stop)
    {
        vector<string> terms;
        unsigned int left_term, right_term;
        left_term = (unsigned int)floor(start * TIME_TERMS_PER_DAY);
        right_term = (unsigned int)ceil(stop * TIME_TERMS_PER_DAY);

        for (unsigned int i = left_term; i < right_term; i++)
            terms.push_back(std::to_string(i));

        return terms;
    }

    double position_angle(const S2Point &a, const S2Point &b)
    {
        // Meeus 1998, Astronomical Algorithms, p116
        S2LatLng aa(a);
        S2LatLng bb(b);
        S1Angle dra = bb.lng() - aa.lng();
        return atan2(sin(dra), cos(aa.lat()) * tan(bb.lat()) - sin(aa.lat()) * cos(dra));
    }

    vector<string> split(string str, const char delimiter)
    {
        int start = 0, end;
        vector<string> parts;
        while ((end = str.find(delimiter, start)) != string::npos)
        {
            parts.push_back(str.substr(start, end - start));
            start = end + 1;
        }
        parts.push_back(str.substr(start)); // remainder of string
        return parts;
    }

    string join(vector<string> s, const char *delimiter)
    {
        if (s.size() == 0)
            return "";

        return std::accumulate(std::next(s.begin()), s.end(), s[0],
                               [delimiter](string a, string b)
                               { return std::move(a) + delimiter + std::move(b); });
    }

    string format_vertices(vector<S2LatLng> vertices)
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

    string format_vertices(vector<S2Point> vertices)
    {
        vector<S2LatLng> ll_vertices;
        for (auto vertex : vertices)
            ll_vertices.push_back(S2LatLng(vertex));
        return format_vertices(ll_vertices);
    }

    string format_vertices(S2LatLngRect fov)
    {
        vector<S2LatLng> vertices;
        for (int i = 0; i < 4; i++)
            vertices.push_back(fov.GetVertex(i));
        return format_vertices(vertices);
    }

    string format_vertices(int num_vertices, double *ra, double *dec)
    {
        vector<S2LatLng> vertices;
        for (int i = 0; i < num_vertices; i++)
            vertices.push_back(S2LatLng::FromDegrees(dec[i], ra[i]));
        return format_vertices(vertices);
    }

    vector<S2Point> makeVertices(string fov)
    {
        vector<S2Point> vertices;
        for (auto coord : split(fov, ','))
        {
            vector<string> values = split(coord, ':');
            if (values.size() < 2)
            {
                std::cerr << "\nfov error: \"" << fov << '"' << std::endl;
                throw std::runtime_error("Not enough vertices");
            }
            try
            {
                S2LatLng ll = S2LatLng::FromDegrees(std::stod(values[1]), std::stod(values[0]));
                vertices.push_back(ll.ToPoint());
            }
            catch (std::invalid_argument const &ex)
            {
                std::cerr << "\nfov error: " << fov << '"' << std::endl;
                throw std::runtime_error("Could not parse fov into vertices");
            }
        }
        return vertices;
    }

    std::unique_ptr<S2Polygon> makePolygon(vector<S2Point> vertices)
    {
        int n;
        n = vertices.size();

        S2Builder::Options builder_options;
        builder_options.set_split_crossing_edges(true);
        S2Builder builder{builder_options};

        S2Polygon polygon;

        s2builderutil::S2PolygonLayer::Options layer_options;
        layer_options.set_edge_type(S2Builder::EdgeType::UNDIRECTED);
        builder.StartLayer(std::make_unique<s2builderutil::S2PolygonLayer>(&polygon, layer_options));

        for (int i = 1; i < n; i++)
            builder.AddEdge(vertices[i - 1], vertices[i]);

        // close the polygon
        builder.AddEdge(vertices[n - 1], vertices[0]);

        S2Error error;
        builder.Build(&error);
        if (!error.ok())
        {
            std::cerr << error.code() << " " << error.text() << std::endl;
            throw std::runtime_error("Polygon build error");
        }

        std::unique_ptr<S2Polygon> result;
        vector<std::unique_ptr<S2Loop>> loops = polygon.Release();
        result = std::make_unique<S2Polygon>(std::move(loops));

        return std::move(result);
    }

    std::unique_ptr<S2Polygon> makePolygon(string fov)
    {
        vector<S2Point> vertices = makeVertices(fov);
        return makePolygon(vertices);
    }

}