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
using std::string;
using std::vector;

namespace sbsearch
{
    namespace testing
    {
        bool is_increasing(const vector<double> &v)
        {
            auto i = std::adjacent_find(v.begin(), v.end(), std::greater<double>());
            return (i == v.end());
        }

    }
}