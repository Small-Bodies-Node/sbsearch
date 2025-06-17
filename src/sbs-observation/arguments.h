#ifndef SBS_OBSERVATION_ARGUMENTS_H_
#define SBS_OBSERVATION_ARGUMENTS_H_

#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include "../cli.h"
#include "../date.h"

using std::optional;
using std::string;
using std::vector;

namespace sbsearch::sbs_observation
{
    // Input file/stream format
    enum FileFormat
    {
        JSON,
        CSV,
        AUTO // determine format from file extension
    };

    // CLI arguments
    struct Arguments : sbsearch::cli::CommonArguments
    {
        string action;
        string file;

        FileFormat file_format;
        int batch_size;
        bool drop_indices;
        bool noop;

        vector<string> sources;
        optional<Date> start_date, stop_date;
    };

    // Read file format from a stream.
    std::istream &operator>>(std::istream &in, FileFormat &file_format);

    // Get arguments from CLI.
    Arguments get_arguments(int argc, char *argv[]);
}

#endif