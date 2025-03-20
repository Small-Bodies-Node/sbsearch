#include "config.h"

#include <cinttypes>
#include <string>
#include <sstream>
#include <vector>

#include "moving_target.h"
#include "util.h"

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

    MovingTarget::MovingTarget(const string &designation, const int moving_target_id, const bool small_body)
    {
        designation_ = designation;
        moving_target_id_ = moving_target_id;
        small_body_ = small_body;
    }

    MovingTarget::MovingTarget(const MovingTarget &other)
    {
        designation_ = other.designation();
        alternate_names_ = set<string>(other.alternate_names());
        moving_target_id_ = other.moving_target_id();
        small_body_ = other.small_body();
    }

    bool MovingTarget::operator==(const MovingTarget &other) const
    {
        return ((moving_target_id_ == other.moving_target_id()) &
                (designation_ == other.designation()) &
                (alternate_names_ == other.alternate_names()) &
                (small_body_ == other.small_body()));
    }

    bool MovingTarget::operator!=(const MovingTarget &other) const
    {
        return !(*this == other);
    }

    std::ostream &operator<<(std::ostream &os, const MovingTarget &target)
    {
        os << to_string(target);
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
            if (designation_ != "")
                alternate_names_.insert(designation_);

            designation_ = name;
        }
        else
            alternate_names_.insert(name);
    }

    string to_string(const MovingTarget &target)
    {
        char s[1024];
        const std::vector<string> alternate_names(target.alternate_names().begin(), target.alternate_names().end());
        sprintf(s, "%s (ID=%" PRId64 "; %ssmall body=%s)",
                target.designation().c_str(),
                target.moving_target_id().value_or(-1),
                alternate_names.size() > 0 ? (join(alternate_names, ", ") + "; ").c_str() : "",
                target.small_body() ? "true" : "false");
        return string(s);
    }
}