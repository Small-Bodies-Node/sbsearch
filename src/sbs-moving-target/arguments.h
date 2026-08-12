#ifndef SBS_MOVING_TARGET_ARGUMENTS_H_
#define SBS_MOVING_TARGET_ARGUMENTS_H_

#include <string>
#include <vector>

#include "config.h"
#include "date.h"
#include "cli.h"

using namespace sbsearch;
using namespace sbsearch::cli;
using std::optional;
using std::string;
using std::vector;

namespace sbsearch::sbs_moving_target
{
    struct Arguments : CommonArguments
    {
        string action;
        string target;
        vector<string> alternate_names;
        bool force_remove;
        optional<Date> start_date, stop_date;
        bool major_body;
    };

    Arguments get_arguments(int argc, char *argv[]);
}

#endif