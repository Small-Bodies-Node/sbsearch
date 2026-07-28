#include "config.h"
#include <string>
#include "ephemeris.h"

namespace sbsearch
{
    int Ephemeris::normalize_index(const int k, const int max) const
    {
        if ((k < -max) || (k >= max))
            throw std::runtime_error("Invalid index " + std::to_string(k) +
                                     " given number of elements: " + std::to_string(max));
        return k + ((k >= 0) ? 0 : max);
    }
}