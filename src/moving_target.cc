#include "config.h"

#include <cinttypes>
#include <string>
#include <sstream>
#include <vector>

#include "moving_target.h"
#include "orbital_elements.h"
#include "util/string.h"

using std::optional;
using std::string;

namespace sbsearch
{
    MovingTarget::MovingTarget(const string &designation, const bool small_body)
    {
        designation_ = designation;
        small_body_ = small_body;
    }

    MovingTarget::MovingTarget(const string &designation, const optional<int64_t> moving_target_id, const bool small_body)
    {
        designation_ = designation;
        moving_target_id_ = moving_target_id;
        small_body_ = small_body;
    }

    MovingTarget::MovingTarget(const string &designation, const OrbitalElements &orbit)
    {
        designation_ = designation;
        orbit_ = orbit;
    }

    MovingTarget::MovingTarget(const MovingTarget &other)
    {
        designation_ = other.designation_;
        alternate_names_ = set<string>(other.alternate_names_);
        moving_target_id_ = other.moving_target_id_;
        small_body_ = other.small_body_;
        orbit_ = other.orbit_;
    }

    bool MovingTarget::operator==(const MovingTarget &other) const
    {
        return ((moving_target_id_ == other.moving_target_id_) &&
                (designation_ == other.designation_) &&
                (alternate_names_ == other.alternate_names_) &&
                (small_body_ == other.small_body_) &&
                (orbit_ == other.orbit_));
    }

    bool MovingTarget::operator!=(const MovingTarget &other) const
    {
        return !(*this == other);
    }

    string MovingTarget::to_string() const
    {
        string alt_names = util::join<string>({alternate_names_.begin(), alternate_names_.end()}, ", ");
        string id = moving_target_id_ ? std::to_string(moving_target_id_.value()) : "null";
        string s = designation_ + " (ID=" + id + "; " +
                   (alt_names.size() == 0 ? "" : alt_names + "; ") +
                   "small body=" + (small_body_ ? "true" : "false") + "; " +
                   "orbit=" + (orbit_ ? "true" : "false") + ")";
        return s;
    }

    std::ostream &operator<<(std::ostream &os, const MovingTarget &target)
    {
        os << target.to_string();
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
            // if already designated, move designation to the alternates
            if (designation_ != "" && designation_ != name)
                alternate_names_.insert(designation_);

            designation_ = name;
        }
        else
            alternate_names_.insert(name);
    }
}