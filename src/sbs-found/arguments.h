#ifndef SBS_SBSFOUND_ARGUMENTS_H_
#define SBS_SBSFOUND_ARGUMENTS_H_

#include <optional>
#include <string>
#include <vector>

#include "config.h"
#include "date.h"
#include "cli.h"

using namespace sbsearch::cli;
using std::optional;
using std::string;
using std::vector;

namespace sbsearch::sbs_found
{
    struct Arguments : CommonArguments
    {
        string action;

        string target;
        bool major_body;
        bool input_file;

        vector<string> sources; // not yet implemented
        Date start_date, stop_date;

        string output_filename;
        OutputFormat output_format = TABLE;
        Date::Format date_format = Date::Format::MJD;

        bool force;
    };

    Arguments get_arguments(int argc, char *argv[]);
}

#endif