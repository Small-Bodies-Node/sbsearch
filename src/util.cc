#include "config.h"
#include "util.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <ostream>
#include <memory>
#include <string>
#include <vector>

#include <s2/s1angle.h>
#include <s2/s2builder.h>
#include <s2/s2builderutil_s2polygon_layer.h>
#include <s2/s2latlng.h>
#include <s2/s2loop.h>
#include <s2/s2point.h>
#include <sqlite3.h>

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
    double position_angle(const S2Point &a, const S2Point &b)
    {
        // Meeus 1998, Astronomical Algorithms, p116
        S2LatLng aa(a);
        S2LatLng bb(b);
        S1Angle dra = bb.lng() - aa.lng();
        return atan2(sin(dra), cos(aa.lat()) * tan(bb.lat()) - sin(aa.lat()) * cos(dra));
    }

    S2LatLng offset_by(const S2LatLng &origin, const S1Angle &position_angle, const S1Angle &distance)
    {
        // Spherical trig based on astropy.coordinates.angle_utilities.offset_by()

        // Spherical trigonometry with:
        // Triangle:
        //   A: north pole (or really-really close to it)
        //   B: point
        //   C: result
        // With angles:
        //   A: change in longitude
        //   B: position angle
        //   C:
        // And sides:
        //   a: distance
        //   b: final co-latitude
        //   c: starting co-latitude

        double cos_a = cos(distance);
        double sin_a = sin(distance);
        double cos_c = sin(origin.lat());
        double sin_c = cos(origin.lat());
        double cos_B = cos(position_angle);
        double sin_B = sin(position_angle);

        double cos_b = cos_c * cos_a + sin_c * sin_a * cos_B;
        double xsin_A = sin_a * sin_B * sin_c;
        double xcos_A = cos_a - cos_b * cos_c;

        S1Angle A = S1Angle::Radians(std::atan2(xsin_A, xcos_A));
        bool small_sin_c = sin_c < 1e-12;
        if (small_sin_c)
        {
            A = S1Angle::Radians(PI_2 + cos_c * (PI_2 - position_angle.radians()));
        }

        S1Angle lon = origin.lng() + A;
        S1Angle lat = S1Angle::Radians(std::asin(cos_b));

        return S2LatLng(lat, lon);
    }

    vector<S2LatLng> ellipse(const int n, const S2LatLng &center, const double &a, const double &b, const double &theta)
    {
        vector<S2LatLng> e;
        const S1Angle th = S1Angle::Radians(theta);

        assert(n >= 4);
        for (int i = 0; i < n; i++)
        {
            const S1Angle phi = S1Angle::Radians(2 * PI * i / n);
            const S1Angle rho = S1Angle::Radians(
                a * b / std::sqrt(std::pow(b * cos(phi), 2) + std::pow(a * sin(phi), 2)));
            e.push_back(offset_by(center, th + phi, rho));
        }

        return e;
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
                std::cerr << "fov error on " << coord << " from " << fov << std::endl;
                throw std::runtime_error("Not enough vertices");
            }
            try
            {
                S2LatLng ll = S2LatLng::FromDegrees(std::stod(values[1]), std::stod(values[0]));
                vertices.push_back(ll.ToPoint());
            }
            catch (std::invalid_argument const &ex)
            {
                std::cerr << "fov error on " << coord << " from " << fov << std::endl;
                throw std::runtime_error("Could not parse fov into vertices");
            }
        }
        return vertices;
    }

    void makePolygon(const vector<S2Point> &vertices, S2Polygon &polygon)
    {
        int n;
        n = vertices.size();

        S2Builder::Options builder_options;
        builder_options.set_split_crossing_edges(true);
        S2Builder builder{builder_options};

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
    }

    void makePolygon(string fov, S2Polygon &polygon)
    {
        makePolygon(makeVertices(fov), polygon);
    }
}