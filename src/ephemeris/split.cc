#include "config.h"

#include "ephemeris/ephemeris.h"
#include "polyline.h"
#include "ephemeris/split.h"

namespace sbsearch::ephemeris
{
    vector<Ephemeris> split(const Ephemeris &eph, const double length, const double time)
    {
        if (eph.num_vertices() <= 1)
            return {};

        vector<Ephemeris> segments;
        segments.reserve(std::ceil(make_polyline(eph).GetLength().degrees() / length));
        double arc = 0, period = 0;
        int start = 0;
        for (int i = 0; i < eph.data().size() - 1; i++)
        {
            Ephemeris segment = eph.segment(i);
            arc += make_polyline(segment).GetLength().degrees();
            period += segment.data(1).mjd.value() - segment.data(0).mjd.value();

            if ((arc >= length) || (period >= time) || (i == eph.num_segments() - 1))
            {
                segments.push_back(eph.slice(start, i + 2));
                arc = 0;
                period = 0;
                start = i + 1;
            }
        }
        return segments;
    }
}