#ifndef SBS_EPHEMERIS_ARGUMENTS_H_
#define SBS_EPHEMERIS_ARGUMENTS_H_

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

namespace sbsearch::sbs_ephemeris
{
    struct Arguments : CommonArguments
    {
        string action;

        string file;
        bool input_file;

        string target;
        bool major_body;
        string observer;
        optional<Date> start_date, stop_date;
        string step_size;

        double interpolate = -1;
        string output_filename;
        OutputFormat output_format = TABLE;
        Date::Format date_format = Date::Format::MJD;

        bool remove_all;
        bool cache = true;
    };

    Arguments get_arguments(int argc, char *argv[]);
}

#endif