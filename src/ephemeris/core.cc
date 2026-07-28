#include "config.h"

#include "ephemeris.h"
#include "moving_target.h"
#include "util/math.h"

namespace sbsearch
{
    Ephemeris::Ephemeris(const MovingTarget &target, const Data &data)
    {
        num_vertices_ = data.size();
        num_segments_ = (num_vertices_ == 0) ? 0 : (num_vertices_ - 1);

        target_ = target;
        data_ = Data(data);

        isValid();
    }

    bool Ephemeris::isValid() const
    {
        if (!util::is_monotonically_increasing(mjd()))
            throw std::runtime_error("mjd must be monotonically increasing.");

        return true;
    }
}