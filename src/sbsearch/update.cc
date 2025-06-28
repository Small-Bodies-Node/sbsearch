#include "observation.h"
#include "sbsdb.h"
#include "sbsearch.h"

using sbsearch::sbsdb::Postgresql;

namespace sbsearch
{
    template <typename SBSDB>
    void SBSearch<SBSDB>::update_observations(Observations &observations)
    {
        index_observations(observations);
        sbsdb::update::observations(&db_, observations);
    }

    template void SBSearch<Postgresql>::update_observations(Observations &);
}