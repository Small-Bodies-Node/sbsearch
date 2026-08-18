#ifndef SBS_SBSMOVING_TARGET_LIST_H_
#define SBS_SBSMOVING_TARGET_LIST_H_

#include <string>
#include <vector>

#include "config.h"
#include "date.h"
#include "cli.h"

using namespace sbsearch;
using std::string;
using std::vector;

namespace sbsearch::sbs_moving_target
{
    template <typename DB>
    void list(SBSearch<DB> &sbs);
}

#endif