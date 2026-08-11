#include <string>

#include "config.h"
#include "ephemeris.h"

using std::to_string;

namespace sbsearch
{
    void check_target_id_(const Ephemeris &a, const Ephemeris &b)
    {
        if (a.target().moving_target_id() != b.target().moving_target_id())
            throw std::runtime_error("Attempted to append an ephemeris with a different object ID: " +
                                     a.target().to_string() + " and " + b.target().to_string());
    }

    void check_time_(const Ephemeris &a, const Ephemeris::Data &b)
    {
        // check that b's time axis follows a
        if (a.num_vertices() > 0 && (a.data().back().mjd > b[0].mjd))
            throw std::runtime_error("Attempting to append ephemeris data with an earlier mjd. Compare" +
                                     to_string(a.data().back().mjd.value()) + " to " +
                                     to_string(b[0].mjd.value()));

        // check that b's time axis is in order
        if (b.size() > 1)
        {
            auto i = std::adjacent_find(b.begin(), b.end(),
                                        [](const auto &left, const auto &right)
                                        { return left.mjd > right.mjd; });
            if (i != b.end())
                throw std::runtime_error("mjd must be monotonically increasing.");
        }
    }

    void Ephemeris::append(const Datum &new_datum)
    {
        check_time_(*this, {new_datum});
        data_.push_back(new_datum);
        num_vertices_ = data_.size();
        num_segments_ = (num_vertices_ == 0) ? 0 : (num_vertices_ - 1);
    }

    void Ephemeris::append(Datum &&new_datum)
    {
        check_time_(*this, {new_datum});
        data_.push_back(std::move(new_datum));
        num_vertices_ = data_.size();
        num_segments_ = (num_vertices_ == 0) ? 0 : (num_vertices_ - 1);
    }

    void Ephemeris::append(const Data &new_data)
    {
        check_time_(*this, new_data);
        std::copy(new_data.begin(), new_data.end(), std::back_inserter(data_));
        num_vertices_ = data_.size();
        num_segments_ = (num_vertices_ == 0) ? 0 : (num_vertices_ - 1);
    }

    void Ephemeris::append(Data &&new_data)
    {
        check_time_(*this, new_data);
        std::move(new_data.begin(), new_data.end(), std::back_inserter(data_));
        num_vertices_ = data_.size();
        num_segments_ = (num_vertices_ == 0) ? 0 : (num_vertices_ - 1);
    }

    void Ephemeris::append(const Ephemeris &eph)
    {
        check_target_id_(*this, eph);
        append(eph.data());
    }

    void Ephemeris::append(Ephemeris &&eph)
    {
        check_target_id_(*this, eph);
        append(std::move(eph.data()));
    }
}