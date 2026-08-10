#include "config.h"

#include "ephemeris.h"
#include "exceptions.h"

namespace sbsearch
{
    Ephemeris Ephemeris::slice(const int start) const
    {
        const int i = normalize_index(start, num_vertices_);
        Data subset(data_.begin() + i, data_.end());
        return Ephemeris(target_, subset);
    }

    Ephemeris Ephemeris::slice(const int start, const int stop) const
    {
        const int i = normalize_index(start, num_vertices_);
        const int j = normalize_index(stop, num_vertices_ + 1);

        if (i > j)
            throw EphemerisError("start cannot be greater than stop.");

        Data subset(data_.begin() + i, data_.begin() + j);
        return Ephemeris(target_, subset);
    }
}