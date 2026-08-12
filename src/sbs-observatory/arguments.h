#ifndef SBS_OBSERVATORY_ARGUMENTS_H_
#define SBS_OBSERVATORY_ARGUMENTS_H_

#include <string>

#include "cli.h"
#include "config.h"
#include "observatory.h"

using namespace sbsearch;
using namespace sbsearch::cli;

using std::string;

namespace sbsearch::sbs_observatory
{
    struct Arguments : CommonArguments
    {
        string action;
        Observatory observatory;
        string output_filename;
    };

    Arguments get_arguments(int argc, char *argv[]);
}

#endif