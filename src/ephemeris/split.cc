#include "config.h"

#include "ephemeris.h"

namespace sbsearch
{
    vector<Ephemeris> Ephemeris::split(double length, double time) const
    {
        if (num_vertices_ <= 1)
            return {};

        vector<Ephemeris> segments;
        segments.reserve(std::ceil(as_polyline().GetLength().degrees() / length));
        double arc = 0, period = 0;
        int start = 0;
        for (int i = 0; i < num_segments_; i++)
        {
            Ephemeris segment_ = segment(i);
            arc += segment_.as_polyline().GetLength().degrees();
            period += segment_.data(1).mjd.value() - segment_.data(0).mjd.value();

            if ((arc >= length) || (period >= time) || (i == num_segments_ - 1))
            {
                segments.push_back(slice(start, i + 2));
                arc = 0;
                period = 0;
                start = i + 1;
            }
        }
        return segments;
    }
}