#ifndef SBS_MOVING_TARGET_H_
#define SBS_MOVING_TARGET_H_

#include <string>
#include <set>
#include "ephemeris.h"

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
        MovingTarget() : designation_("") {}
        // define a moving target with its primary designation
        MovingTarget(const string &designation) : designation_(designation) {}
        // copy
        MovingTarget(const MovingTarget &other) : designation_(other.designation())
        {
            alternate_names_ = set<string>(other.alternate_names());
            object_id_ = other.object_id();
        }

        // get/set primary designation
        const string &designation() const { return designation_; };
        void designation(const string &designation) { designation_ = designation; };

        // get all alternate names
        const set<string> &alternate_names() const { return alternate_names_; }

        // add a name, default is to add an alternate
        void add_name(const string &name, const bool primary = false)
        {
            if (primary)
                designation_ = name;
            else
                alternate_names_.insert(name);
        }

        template <typename ForwardIterator>
        void add_names(const ForwardIterator &begin, const ForwardIterator &end)
        {
            alternate_names_.insert(begin, end);
        }

        // get/set database object ID
        const int &object_id() const { return object_id_; };
        void object_id(const int id) { object_id_ = id; };

    private:
        string designation_;
        set<string> alternate_names_;
        int object_id_ = UNDEF_OBJECT_ID;
    };

}

#endif // SBS_MOVING_TARGET_H_