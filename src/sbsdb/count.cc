#include <string>

#include "count.h"
#include "postgresql.h"

using std::string;

namespace sbsearch::sbsdb::count
{
    template <typename DB>
    int observations(DB &db, const double mjd_start, const double mjd_stop)
    {
        return db.template get_one<int>(
            "SELECT COUNT(*) FROM observations WHERE mjd_start >= $1 AND mjd_stop <= $2",
            mjd_start,
            mjd_stop);
    };

    template <typename DB>
    int observations(DB &db, const string &source, const double mjd_start, const double mjd_stop)
    {
        return db.template get_one<int>(
            "SELECT COUNT(*) FROM observations WHERE source = $1 AND mjd_start >= $2 AND mjd_stop <= $3",
            source,
            mjd_start,
            mjd_stop);
    };

    template int observations(Postgresql &, const double, const double);
    template int observations(Postgresql &, const string &, const double, const double);
}