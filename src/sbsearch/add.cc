#include <string>

#include "ephemeris.h"
#include "exceptions.h"
#include "observation.h"
#include "sbsdb.h"
#include "sbsearch.h"

using sbsearch::sbsdb::Postgresql;

namespace sbsearch
{
    template <typename SBSDB>
    void SBSearch<SBSDB>::add_ephemeris(Ephemeris &eph)
    {
        if (sbsdb::count::ephemeris(&db_, eph.target(), eph.data(0).mjd.value(), eph.data(-1).mjd.value()) != 0)
            throw EphemerisError("data already present in database for target and date range: " +
                                 eph.target().to_string() + ", " +
                                 std::to_string(eph.data(0).mjd.value()) + ", " +
                                 std::to_string(eph.data(-1).mjd.value()));

        sbsdb::add::ephemeris(&db_, eph);
    }

    template <typename SBSDB>
    void SBSearch<SBSDB>::add_observations(Observations &observations)
    {
        index_observations(observations);
        sbsdb::add::observations(&db_, observations);
    }

    template void SBSearch<Postgresql>::add_ephemeris(Ephemeris &);
    template void SBSearch<Postgresql>::add_observations(Observations &);
}