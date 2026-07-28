#ifndef SBS_QUERY_ARGUMENTS_H_
#define SBS_QUERY_ARGUMENTS_H_

#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include "cli.h"
#include "date.h"
#include "intersection.h"

using namespace sbsearch;
using namespace sbsearch::cli;
using std::optional;
using std::string;
using std::vector;

namespace sbsearch::sbs_query
{
    // CLI arguments

    struct Arguments : CommonArguments
    {
        string target;
        bool fixed_target;
        bool input_file;

        vector<string> sources;
        bool no_parallax;
        bool use_uncertainty;
        double padding = 0;
        double arc_length = 10;
        double time_period = 90;
        bool approximate;
        bool save;
        string info_file;
        string output_file;
        OutputFormat output_format;
        DateFormat date_format = DateFormat::MJD;
        bool show_fov = false;

        IntersectionType intersection_type = IntersectsArea;

        string eph_file;
        string orbit_file;
        bool horizons;
        bool major_body;

        string observer;
        optional<Date> start_date, stop_date;
        string time_step;

        bool no_cache;

        bool parallax() const { return !no_parallax; };
        bool cache() const { return !no_cache; };
    };

    // Read file format from a stream.
    std::istream &operator>>(std::istream &in, OutputFormat &output_format);

    // Get arguments from CLI.
    Arguments get_arguments(int argc, char *argv[]);
}

#endif