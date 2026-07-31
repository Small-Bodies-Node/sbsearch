#include "config.h"

#include "ephemeris.h"

namespace sbsearch
{
    void Ephemeris::append(const Datum &new_datum)
    {
        // check that the new point's time is OK
        if (num_vertices_ > 0 && data_.back().mjd > new_datum.mjd)
            throw std::runtime_error("Attempting to append an ephemeris point with an earlier mjd.");
        data_.push_back(new_datum);
        num_vertices_ = data_.size();
        num_segments_ = (num_vertices_ == 0) ? 0 : (num_vertices_ - 1);
    }

    void Ephemeris::append(const Data &new_data)
    {
        // check that new_data's time axis is OK
        if (num_vertices_ > 0 && (data_.back().mjd > new_data[0].mjd))
            throw std::runtime_error("Attempting to append ephemeris data with an earlier mjd.");

        if (new_data.size() > 1)
        {
            auto i = std::adjacent_find(new_data.begin(), new_data.end(),
                                        [](const Datum &a, const Datum &b)
                                        { return a.mjd > b.mjd; });
            if (i != new_data.end())
                throw std::runtime_error("mjd must be monotonically increasing.");
        }

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