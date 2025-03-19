
#include <cinttypes>

#include "observation.h"
#include "sbsdb/get.h"
#include "sbsdb/postgresql.h"

namespace sbsearch::sbsdb
{
    Observation get(Postgresql db, int64_t id)
    {
        return Observation();
    };
}