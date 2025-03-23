#ifndef SBSDB_COUNT_H_
#define SBSDB_COUNT_H_

#include <string>

using std::string;

namespace sbsearch::sbsdb::count
{
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
    int observations(DB &db, const double mjd_start, const double mjd_stop);

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
    int observations(DB &db, const string &source, const double mjd_start, const double mjd_stop);
}

#endif // SBSDB_COUNT_H_
