#include "sbsearch_testing.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include <s2/s2point.h>
#include <s2/s2latlng.h>

using std::cerr;
using std::cout;
using std::endl;
using std::vector;

namespace sbsearch
{
    namespace testing
    {
        template <typename T>
        void print_vector(const char *comment, const vector<T> &v, std::ostream &stream)
        {
            stream << comment;
            for (const auto &n : v)
                stream << n << ' ';
            stream << '\n';
        }
        void print_vector(const char *comment, const vector<double> &v, std::ostream &stream = std::cout)
        {
            return print_vector<double>(comment, v, stream);
        }
        void print_vector(const char *comment, const vector<S2Point> &v, std::ostream &stream = std::cout)
        {
            return print_vector<S2Point>(comment, v, stream);
        }
        void print_vector(const char *comment, const vector<S2LatLng> &v, std::ostream &stream = std::cout)
        {
            return print_vector<S2LatLng>(comment, v, stream);
        }

        bool is_increasing(const vector<double> &v)
        {
            auto i = std::adjacent_find(v.begin(), v.end(), std::greater<double>());
            return (i == v.end());
        }
    }
}