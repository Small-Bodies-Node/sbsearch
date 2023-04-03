#include "config.h"

#include <string>
#include <sstream>

#include "moving_target.h"
#include "util.h"

using std::string;

namespace sbsearch
{
    MovingTarget::MovingTarget(const string &designation)
    {
        designation_ = designation;
    }

    MovingTarget::MovingTarget(const string &designation, const int moving_target_id)
    {
        designation_ = designation;
        moving_target_id_ = moving_target_id;
    }

    MovingTarget::MovingTarget(const MovingTarget &other)
    {
        designation_ = other.designation();
        alternate_names_ = set<string>(other.alternate_names());
        moving_target_id_ = other.moving_target_id();
    }

    bool MovingTarget::operator==(const MovingTarget &other) const
    {
        return ((moving_target_id_ == other.moving_target_id()) &
                (designation_ == other.designation()) &
                (alternate_names_ == other.alternate_names()));
    }

    bool MovingTarget::operator!=(const MovingTarget &other) const
    {
        return !(*this == other);
    }

    std::ostream &operator<<(std::ostream &os, const MovingTarget &target)
    {
        os << target.designation() << " (ID=" << target.moving_target_id();
        auto names = target.alternate_names();
        if (names.size() > 0)
            os << join(vector<string>(names.begin(), names.end()), ", ");
        os << ")";
        return os;
    }

    void MovingTarget::designation(const string &designation)
    {
        // avoid having one name in two places
        alternate_names_.erase(designation);
        designation_ = designation;
    };

    void MovingTarget::add_name(const string &name, const bool primary)
    {
        if (primary)
        {
            alternate_names_.insert(designation_);
            designation_ = name;
        }
        else
            alternate_names_.insert(name);
    }

    string to_string(const MovingTarget &target)
    {
        std::ostringstream sstr;
        sstr << target;
        return sstr.str();
    }
}