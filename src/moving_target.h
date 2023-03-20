#ifndef SBS_MOVING_TARGET_H_
#define SBS_MOVING_TARGET_H_

#include <string>
#include <set>

#define UNDEF_OBJECT_ID -1

using std::set;
using std::string;

namespace sbsearch
{
    // designation == primary identifier, e.g., the name one would use with an
    //   ephemeris service
    // alternate_names == any unique name
    class MovingTarget
    {
    public:
        MovingTarget() {}
        // define a moving target with its primary designation
        MovingTarget(const string &designation);
        // primary designation and object_id
        MovingTarget(const string &designation, const int object_id);
        // copy
        MovingTarget(const MovingTarget &other);

        // strict comparison, must match designation, object_id, and
        // alternate_names
        bool operator==(const MovingTarget &other) const;
        bool operator!=(const MovingTarget &other) const;

        // get primary designation
        inline const string &designation() const { return designation_; };

        // Set primary designation.
        //
        // The previous designation is discarded.  Use add_name to preserve it.
        //
        // If the new name was an alternate name, it is removed from the
        // alternate_name set.
        void designation(const string &designation);

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

        // get/set database object ID
        inline const int &object_id() const { return object_id_; };
        inline void object_id(const int id) { object_id_ = id; };

    private:
        string designation_ = "";
        set<string> alternate_names_;
        int object_id_ = UNDEF_OBJECT_ID;
    };
}

#endif // SBS_MOVING_TARGET_H_