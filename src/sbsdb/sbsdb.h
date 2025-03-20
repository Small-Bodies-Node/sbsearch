#ifndef SBSDB_SBSDB_H_
#define SBSDB_SBSDB_H_

#include <cinttypes>
#include <string>

using std::string;

namespace sbsearch::sbsdb
{
    /**
     * @brief A row in the moving_targets database table.
     *
     */
    struct MovingTargetsModel
    {
        int64_t moving_targets_row_id;
        int64_t moving_target_id;
        string name;
        bool small_body;
        bool primary_id;
    };
}

#endif // SBSDB_SBSDB_H_
