#include "sbsearch_testing.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include <s2/s2latlng.h>
#include <s2/s2loop.h>
#include <s2/s2point.h>
#include <s2/s2point_span.h>

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace sbsearch
{
    namespace testing
    {
        void print_loop(const char *comment, S2Loop *loop, const char *separator, std::ostream &stream)
        {
            S2PointSpan points = loop->vertices_span();
            vector<S2Point> vertices(points.size());
            std::copy(points.begin(), points.end(), vertices.begin());
            print_vertices(comment, vertices, separator, stream);
        }
        void print_vertices(const char *comment, vector<S2Point> vertices, const char *separator, std::ostream &stream)
        {
            vector<S2LatLng> coords(vertices.size());
            std::transform(vertices.begin(), vertices.end(), coords.begin(),
                           [](const S2Point &p)
                           { return S2LatLng(p); });
            print_vector(comment, coords, separator, stream);
        }
    }
}