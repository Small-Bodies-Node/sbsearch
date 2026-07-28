#ifndef SBS_OBSERVATION_ARGUMENTS_H_
#define SBS_OBSERVATION_ARGUMENTS_H_

#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include "cli.h"
#include "date.h"

using namespace sbsearch::cli;
using std::optional;
using std::string;
using std::vector;

namespace sbsearch::sbs_observation
{
    // CLI arguments
    struct Arguments : CommonArguments
    {
        string action;
        string file;

        OutputFormat file_format;
        int batch_size;
        bool drop_indices;
        bool noop;

        vector<string> sources;
        optional<Date> start_date, stop_date;
    };

    // Get arguments from CLI.
    Arguments get_arguments(int argc, char *argv[]);
}

#endif