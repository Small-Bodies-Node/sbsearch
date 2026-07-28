#include "config.h"

#include "ephemeris.h"

namespace sbsearch
{
    int Ephemeris::num_vertices() const
    {
        return num_vertices_;
    }

    S2Point Ephemeris::vertex(const int k) const
    {
        return data(k).as_s2point();
    }

    vector<S2Point> Ephemeris::vertices() const
    {
        vector<S2Point> result(num_vertices_);
        for (int i = 0; i < num_vertices_; i++)
            result[i] = vertex(i);
        return result;
    }

    S2Polyline Ephemeris::as_polyline() const
    {
        return S2Polyline(vertices());
    }
}