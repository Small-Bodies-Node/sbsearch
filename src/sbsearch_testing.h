#ifndef SBSEARCH_TESTING_H_
#define SBSEARCH_TESTING_H_

#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <vector>
#include <s2/s2loop.h>

using std::vector;

namespace sbsearch
{
    namespace testing
    {
        template <typename T>
        void print_vector(const char *comment, const vector<T> &sequence, const char *separator = " ", std::ostream &stream = std::cout);
        void print_loop(const char *comment, S2Loop *loop, const char *separator = " ", std::ostream &stream = std::cout);
        void print_vertices(const char *comment, vector<S2Point> vertices, const char *separator = " ", std::ostream &stream = std::cout);

        bool is_increasing(const vector<double> &v);

        // define templates
        template <typename T>
        void print_vector(const char *comment, const vector<T> &v, const char *separator, std::ostream &stream)
        {
            stream << comment;
            for (const auto &n : v)
                stream << n << separator;
            stream << '\n';
        }
    }
}

#endif
