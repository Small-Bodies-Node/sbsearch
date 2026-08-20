#include <iostream>
#include <stdexcept>
#include <string>
#include <boost/program_options.hpp>

#include "config.h"
#include "date.h"
#include "cli.h"
#include "sbs-ephemeris/arguments.h"

using namespace sbsearch;
using namespace sbsearch::cli;
using std::cerr;
using std::cout;
using std::endl;
using std::string;

namespace sbsearch::sbs_ephemeris
{
    Arguments get_arguments(int argc, char *argv[])
    {
        using namespace boost::program_options;

        Arguments args;

        positional_options_description positional;
        positional.add("action", 1).add("target", 1);

        options_description hidden("Hidden options");
        hidden.add_options()("action", value<string>(&args.action), "ephemeris action");
        hidden.add_options()("target", value<string>(&args.target), "target");

        options_description target_options("Target options");
        target_options.add_options()(
            "input,i", bool_switch(&args.input_file), "read target names from an input file");

        options_description date_range("Options for date ranges");
        date_range.add_options()(
            "start", value<optional<Date>>(&args.start_date),
            "start date for ephemeris data [YYYY-MM-DD or MJD]")(
            "stop,end", value<optional<Date>>(&args.stop_date),
            "stop date for ephemeris data [YYYY-MM-DD or MJD]")(
            "step", value<string>(&args.step_size)->default_value("auto"), "time step for Horizons ephemeris generation: length and time unit, \"auto\" for a variable step size based on distance, or \"VAR X\" where X is an angular distance in arcsec (use with caution)");

        options_description add_options("Options for add action");
        add_options.add_options()(
            "file", value<string>(&args.file), "read ephemeris from this file (JSON or Horizons format)")(
            "no-cache", value<bool>(&args.cache)->implicit_value(false), "do not allow cached Horizons results")(
            "major-body", bool_switch(&args.major_body), "moving target is a major body");

        options_description clean_cache_options("Options for clean-cache action");
        clean_cache_options.add_options()(
            "max-age", value<double>(&args.max_age)->default_value(3), "remove cached data older than max-age days");

        options_description get_options("Options for get action");
        get_options.add_options()(
            "output,o", value<string>(&args.output_filename), "save ephemeris to this file")(
            "observer", value<string>(&args.observer)->default_value("500@399"), "observer location")(
            "no-cache", value<bool>(&args.cache)->implicit_value(false), "do not allow cached Horizons results")(
            "format,f", value<OutputFormat>(&args.output_format), "output file format: table (default) or json")(
            "date", value<Date::Format>(&args.date_format), "date format: mjd (default) or calendar");

        options_description list_options("Options for list action");
        list_options.add_options()(
            "interpolate", value<double>(&args.interpolate), "interpolate the database ephemeris to this time step in days")(
            "observer", value<string>(&args.observer)->default_value("500@399"), "offset ephemeris for this observer location")(
            "output,o", value<string>(&args.output_filename), "save ephemeris to this file")(
            "format,f", value<OutputFormat>(&args.output_format), "output file format: table (default) or json")(
            "date", value<Date::Format>(&args.date_format), "date format: mjd (default) or calendar");

        options_description remove_options("Options for remove action");
        remove_options.add_options()(
            "all", bool_switch(&args.remove_all), "remove all ephemeris data");

        options_description general = get_common_options((CommonArguments *)&args);

        options_description visible("");
        visible
            .add(target_options)
            .add(date_range)
            .add(add_options)
            .add(clean_cache_options)
            .add(list_options)
            .add(remove_options)
            .add(general);

        options_description all("");
        all.add(visible).add(hidden);

        options_description add_action("");
        add_action.add(add_options).add(date_range).add(general);

        options_description clean_cache_action("");
        clean_cache_action.add(clean_cache_options).add(general);

        options_description get_action("");
        get_action.add(date_range).add(get_options).add(general);

        options_description list_action("");
        list_action.add(date_range).add(list_options).add(general);

        options_description remove_action("");
        remove_action.add(remove_options).add(date_range).add(general);

        variables_map vm;
        boost::program_options::store(command_line_parser(argc, argv).options(all).positional(positional).run(), vm);
        boost::program_options::notify(vm);

        if (vm.count("version"))
        {
            cout << "SBSearch " << SBSEARCH_VERSION << "\n";
            exit(0);
        }

        if (vm.count("help") || !vm.count("action"))
        {
            // help for a specific action?
            if (args.action == "add")
            {
                cout << "Usage: sbs-ephemeris add <target> [options...]\n"
                     << "Add ephemeris data for a target to the database.\n\n"
                     << "<target> is the ephemeris target name / designation\n"
                     << "or, with -i, an input file listing multiple targets\n"
                     << add_action << "\n";
            }
            else if (args.action == "clean-cache")
            {
                cout << "Usage: sbs-ephemeris clean-cache [options...]\n"
                     << "Clean the ephemeris data cache.\n\n"
                     << clean_cache_action << "\n";
            }
            else if (args.action == "get")
            {
                cout << "Usage: sbs-ephemeris get <target> [options...]\n"
                     << "Get ephemeris data from Horizons.\n\n"
                     << get_action << "\n";
            }
            else if (args.action == "list")
            {
                cout << "Usage: sbs-ephemeris list <target> [options...]\n"
                     << "List ephemeris data for a target in the database.\n\n"
                     << list_action << "\n";
            }
            else if (args.action == "remove")
            {
                cout << "Usage: sbs-ephemeris remove <target> [options...]\n"
                     << "Remove ephemeris data for a target from the database.\n\n"
                     << "<target> is the ephemeris target name / designation\n"
                     << remove_action << "\n";
            }
            else
            {
                cout << "Usage: sbs-ephemeris <action> <target> [options...]\n\n"
                     << "Manage sbsearch ephemeris data.\n\n"
                     << "<action> is one of {add, clean-cache, get, list, remove}\n"
                     << "<target> is the ephemeris target name / designation\n"
                     << visible << "\n";
            }

            if ((args.action == "add") || (args.action.empty()))
            {
                cout << "Horizons ephemeris files require the CSV format, angles formatted in degrees,\n"
                     << "dates as Julian days, range in au.  Use the ICRF reference frame and\n"
                     << "observer quantities 1, 9, 19, 20, 23, 24, 27, 37, 41, and 47.  Extra\n"
                     << "precision is optional.\n";
            }

            if (!vm.count("action"))
                cout << "\naction is a required argument\n";

            exit(0);
        }

        if (vm.count("action") && (args.action != "clean-cache") && !vm.count("target"))
            cout << "\target is a required argument\n";

        validate_common_options(vm);
        action_is(args.action, {"add", "clean-cache", "get", "list", "remove"});
        conflicting_options(vm, "file", "horizons");
        conflicting_options(vm, "file", "input");
        option_dependency(vm, "horizons", "start");
        option_dependency(vm, "start", "stop");

        if ((args.action == "add") && (args.observer != "500@399"))
            throw std::logic_error("add action and --observer are not compatible");

        if ((args.action == "get") && (args.input_file))
            throw std::logic_error("get action and --input-file are not compatible");

        if ((args.action == "list") && (args.input_file))
            throw std::logic_error("list action and --input-file are not compatible");

        if ((args.action == "remove") && (!args.remove_all) && (!vm.count("start")))
            throw std::logic_error("remove action requires a date range or --all");

        if (args.start_date && args.stop_date && (args.start_date.value() >= args.stop_date.value()))
            throw std::logic_error("start_date must be before stop_date");

        return args;
    }

}