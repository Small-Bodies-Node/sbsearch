#ifndef SBS_MOVING_TARGET_H_
#define SBS_MOVING_TARGET_H_

#include <optional>
#include <string>
#include <string_view>
#include <set>

#include "orbital_elements.h"

using std::optional;
using std::set;
using std::string;
using std::string_view;

namespace sbsearch
{
    // designation == primary identifier, e.g., the name one would use with an
    //   ephemeris service
    // alternate_names == any unique name
    class MovingTarget
    {
    public:
        MovingTarget() {}

        // define a moving target with its primary designation and small body flag
        MovingTarget(string_view designation, const bool small_body = true);

        // primary designation, moving_target_id, and small body flag
        MovingTarget(string_view designation,
                     const optional<int64_t> moving_target_id,
                     const bool small_body);

        // name and orbit
        MovingTarget(string_view designation, const OrbitalElements &orbit);

        // copy
        MovingTarget(const MovingTarget &other);

        // strict comparison, must match designation, moving_target_id,
        // small_body flag, alternate_names, and orbit
        bool operator==(const MovingTarget &other) const;
        bool operator!=(const MovingTarget &other) const;

        // Describe the target as a string.
        string to_string() const;

        friend std::ostream &operator<<(std::ostream &os, const MovingTarget &target);

        // get primary designation
        inline string_view designation() const { return designation_; };

        // Set primary designation.
        //
        // The previous designation is discarded.  Use add_name to preserve it.
        //
        // If the new name was an alternate name, it is removed from the
        // alternate_name set.
        void designation(const string &designation);

        // get/set small_body flag
        inline const bool &small_body() const { return small_body_; }
        void small_body(const bool flag) { small_body_ = flag; }

        // get all alternate names
        inline const set<string> &alternate_names() const { return alternate_names_; }

        // Add a name.
        //
        // The default is to add an alternate.  If primary is true, then move
        // old designation to an alternate and update the designation with the
        // new name.
        void add_name(const string &name, const bool primary = false);

        template <typename ForwardIterator>
        void add_names(const ForwardIterator &begin, const ForwardIterator &end)
        {
            alternate_names_.insert(begin, end);
        }

        // get/set database moving target ID
        inline const optional<int64_t> &moving_target_id() const { return moving_target_id_; };
        inline void moving_target_id(const optional<int64_t> id) { moving_target_id_ = id; };

        // get/set orbit
        inline const optional<OrbitalElements> &orbit() const { return orbit_; };
        inline void orbit(const optional<OrbitalElements> orbit) { orbit_ = orbit; };

        // A row in the moving_targets database table.
        struct DBModel
        {
            int64_t moving_targets_row_id;
            int64_t moving_target_id;
            string name;
            bool small_body;
            bool primary_id;
        };

    private:
        string designation_ = "";
        set<string> alternate_names_;
        optional<int64_t> moving_target_id_;
        bool small_body_ = true;
        optional<OrbitalElements> orbit_ = std::nullopt;
    };
}

#endif // SBS_MOVING_TARGET_H_