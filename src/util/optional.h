#ifndef SBSEARCH_UTIL_OPTIONAL_H_
#define SBSEARCH_UTIL_OPTIONAL_H_

#include <optional>
#include <vector>

using std::optional;
using std::vector;

namespace sbsearch::util
{

    // Realize a vector of optional values to their values, or else throw
    // std::bad_optional_access
    template <typename T>
    vector<T> optionals_to_values(const vector<optional<T>> vo)
    {
        vector<T> v;
        v.reserve(vo.size());
        std::transform(vo.cbegin(), vo.cend(), std::back_inserter(v),
                       [](const auto &value)
                       { return value.value(); });
        return v;
    };
}
#endif
