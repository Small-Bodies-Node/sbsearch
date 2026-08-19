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
                       [](const auto &item)
                       { return item.value(); });
        return v;
    };

    // Realize a vector of optional values to their values, or else use a default.
    template <typename T>
    vector<T> optionals_to_values(const vector<optional<T>> vo, T default_value)
    {
        vector<T> v;
        v.reserve(vo.size());
        std::transform(vo.cbegin(), vo.cend(), std::back_inserter(v),
                       [&default_value](const auto &item)
                       { return item.value_or(default_value); });
        return v;
    };
}
#endif
