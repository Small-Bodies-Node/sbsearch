#include "config.h"

#include "ephemeris.h"

namespace sbsearch
{
    void Ephemeris::append(const Data &new_data)
    {
        if (num_vertices_ != 0)
            if (data(-1).mjd > new_data[0].mjd)
                throw std::runtime_error("Attempting to append an ephemeris with an earlier mjd.");

        // check that new_data's time axis is OK
        auto i = std::adjacent_find(new_data.begin(), new_data.end(),
                                    [](const Datum &a, const Datum &b)
                                    { return a.mjd > b.mjd; });
        if (i != new_data.end())
            throw std::runtime_error("mjd must be monotonically increasing.");

        data_.insert(data_.end(), new_data.begin(), new_data.end());
        num_vertices_ = data_.size();
        num_segments_ = (num_vertices_ == 0) ? 0 : (num_vertices_ - 1);
    }

    void Ephemeris::append(const Ephemeris &eph)
    {
        if (eph.target().moving_target_id() != target_.moving_target_id())
            throw std::runtime_error("Attempted to append an ephemeris with a different object ID.");

        append(eph.data());
    }
}