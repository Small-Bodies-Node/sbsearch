#include <string>

#include "moving_target.h"
#include "sbsdb/count.h"
#include "sbsdb/postgresql.h"

using std::string;

namespace sbsearch::sbsdb::count
{
    template <typename DB>
    int64_t ephemeris(DB *db, const MovingTarget &target, const double mjd_start, const double mjd_stop)
    {
        if (!target.moving_target_id())
            return 0;

        return db->template get_one<int64_t>(
            "SELECT COUNT(*) FROM ephemerides WHERE moving_target_id=$1 AND mjd >= $2 AND mjd <= $3",
            target.moving_target_id().value(),
            mjd_start,
            mjd_stop);
    }

    template <typename DB>
    int64_t observations(DB *db, const double mjd_start, const double mjd_stop)
    {
        return db->template get_one<int64_t>(
            "SELECT COUNT(*) FROM observations WHERE mjd_start >= $1 AND mjd_stop <= $2",
            mjd_start,
            mjd_stop);
    };

    template <typename DB>
    int64_t observations(DB *db, string_view source, const double mjd_start, const double mjd_stop)
    {
        return db->template get_one<int64_t>(
            "SELECT COUNT(*) FROM observations WHERE source = $1 AND mjd_start >= $2 AND mjd_stop <= $3",
            source,
            mjd_start,
            mjd_stop);
    };

    template int64_t ephemeris(Postgresql *, const MovingTarget &, const double, const double);
    template int64_t observations(Postgresql *, const double, const double);
    template int64_t observations(Postgresql *, string_view, const double, const double);
}