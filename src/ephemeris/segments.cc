#include "config.h"

#include "ephemeris.h"

namespace sbsearch
{
    int Ephemeris::num_segments() const
    {
        return num_segments_;
    }

    Ephemeris Ephemeris::segment(const int k) const
    {
        const int i = normalize_index(k, num_segments_);
        return Ephemeris(target_, {data_[i], data_[i + 1]});
    }

    vector<Ephemeris> Ephemeris::segments() const
    {
        vector<Ephemeris> eph(num_segments_);
        for (int i = 0; i < num_segments_; i++)
            eph[i] = segment(i);
        return eph;
    }
}