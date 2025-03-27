#ifndef SBSDB_COUNT_H_
#define SBSDB_COUNT_H_

#include <string>

#include "../moving_target.h"

using std::string;

namespace sbsearch::sbsdb::count
{
    /**
     * @brief Count ephemeris vertices for a moving target within a time range.
     *
     * @param db The sbsearch database.
     *
     * @param target The moving target (ID).
     *
     * @param mjd_start The start modified Julian date.
     *
     * @param mjd_stop The end modified Julian date.
     *
     * @return int
     */
    template <typename DB>
    int64_t ephemeris(DB &db, const MovingTarget &target, const double mjd_start = 0, const double mjd_stop = 100000);

    /**
     * @brief Count observations within a time range.
     *
     * @param db The sbsearch database.
     *
     * @param mjd_start The start modified Julian date.
     *
     * @param mjd_stop The end modified Julian date.
     *
     * @return int
     */
    template <typename DB>
    int64_t observations(DB &db, const double mjd_start, const double mjd_stop);

    /**
     * @brief Count observations within a time range for a single source.
     *
     * @param db The sbsearch database.
     *
     * @param source The name of the data source.
     *
     * @param mjd_start The start modified Julian date.
     *
     * @param mjd_stop The end modified Julian date.
     *
     * @return int
     */
    template <typename DB>
    int64_t observations(DB &db, const string &source, const double mjd_start, const double mjd_stop);
}

#endif // SBSDB_COUNT_H_
