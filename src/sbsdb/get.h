#include <cinttypes>

namespace sbsearch::sbsdb
{
    template <typename DB, typename T>
    T get(DB db, int64_t id);
};
