#ifndef SBSEARCH_TESTING_H_
#define SBSEARCH_TESTING_H_

#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;
using std::vector;

namespace sbsearch
{
    namespace testing
    {
        template <typename T>
        void print_vector(const char *comment, const vector<T> &sequence, std::ostream &stream = std::cout);

        bool is_increasing(const vector<double> &v);

        bool is_less_than_zero(double x);
    }
}

#endif