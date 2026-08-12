#include <optional>
#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "config.h"
#include "date.h"
#include "./arguments.h"

using namespace sbsearch;
using namespace sbsearch::cli;
using std::cerr;
using std::cout;
using std::endl;
using std::optional;
using std::string;
using std::vector;

namespace sbsearch::sbs_query
{
    Arguments get_arguments(int argc, char *argv[])
    {
        using namespace boost::program_options;

        Arguments args;

        positional_options_description positional;
        positional.add("target", 1);

        options_description hidden("Hidden options");
        hidden.add_options()("target", value<string>(&args.target), "target");

        options_description common_options("Moving / fixed target common options");
        common_options.add_options()(
            "input,i", bool_switch(&args.input_file), "read target names from an input file")(
            "fixed", bool_switch(&args.fixed_target), "indicates <target> is an RA, Dec pair in degrees, e.g., \"123.45 67.890\"")(
            "source,s", value<vector<string>>(&args.sources), "only search this source data set, may be specified multiple times")(
            "padding,p", value<double>(&args.padding), "areal search around query, in arcminutes")(
            "arc-length,arc", value<double>(&args.arc_length), "maximum arc length for ephemeris splitting, degrees")(
            "time-period", value<double>(&args.time_period), "maximum time period for ephemeris splitting, days")(
            "approximate,approx", bool_switch(&args.approximate), "return approximate results")(
            "output,o", value<string>(&args.output_file), "save the results to this file")(
            "format,f", value<OutputFormat>(&args.output_format)->default_value(OutputFormat::AUTO), "output file format: table or json; default is based on the suffix")(
            "date", value<Date::Format>(&args.date_format)->default_value(Date::Format::MJD), "table date format: mjd (default) or calendar")(
            "show-fov", bool_switch(&args.show_fov), "show fields of view in output table")(
            "info", value<string>(&args.info_file), "save query information to this file, JSON format");

        options_description fixed_target_options("Fixed target options");
        fixed_target_options.add_options()(
            "intersection-type", value<IntersectionType>(&args.intersection_type), "areal search type: ContainsArea, IntersectsArea (default), or ContainedByArea");

        options_description moving_target_options("Moving target options");
        moving_target_options.add_options()(
            "all", bool_switch(&args.all_moving_targets), "search all moving targets")(
            "major-body", bool_switch(&args.major_body), "moving target is a major body")(
            "format-help", "display help on file formats and exit")(
            "eph-file", value<string>(&args.eph_file), "read ephemeris from this file (JSON or Horizons format)")(
            "orbit-file", value<string>(&args.orbit_file), "read orbital elements from this file (JSON format)")(
            "horizons", bool_switch(&args.horizons), "generate ephemeris with JPL/Horizons")(
            "start", value<optional<Date>>(&args.start_date), "start date for query [YYYY-MM-DD or MJD]")(
            "stop,end", value<optional<Date>>(&args.stop_date), "stop date for query [YYYY-MM-DD or MJD]")(
            "step", value<string>(&args.step_size)->default_value("auto"), "time step for Horizons ephemeris generation: length and time unit, \"auto\" for a variable step size based on distance, or \"VAR X\" where X is an angular distance in arcsec (use with caution)")(
            "use-uncertainty,u", bool_switch(&args.use_uncertainty), "areal search around ephemeris position using the ephemeris uncertainty")(
            "no-cache", value<bool>(&args.cache)->implicit_value(false), "do not use a file cache for Horizons results")(
            "no-parallax", bool_switch(&args.no_parallax), "do not account for moving target parallax between observatory and the Earth's center")(
            "save", bool_switch(&args.save), "save to/update the found object database");

        options_description general = get_common_options((CommonArguments *)&args);

        options_description visible("");
        visible.add(common_options)
            .add(fixed_target_options)
            .add(moving_target_options)
            .add(general);

        options_description all("");
        all.add(visible).add(hidden);

        variables_map vm;
        boost::program_options::store(
            command_line_parser(argc, argv)
                .options(all)
                .positional(positional)
                .run(),
            vm);
        boost::program_options::notify(vm);

        if (vm.count("version"))
        {
            cout << "SBSearch " << SBSEARCH_VERSION << "\n";
            exit(0);
        }

        if (vm.count("help") || (!vm.count("target") && !vm.count("all")))
        {
            cout << "Usage: sbs-query <target> [options...]\n\n"
                 << "Find observations of a moving or fixed target.\n\n"
                 << "<target> is the ephemeris target name / designation or,\n"
                 << "with --fixed-target, an RA, Dec pair.  Use --input to\n"
                 << "indicate that <target> is a file listing multiple\n"
                 << "targets\n"
                 << visible << "\n";

            if (!vm.count("target"))
                cout << "\ntarget is a required argument\n";

            exit(0);
        }

        validate_common_options(vm);
        conflicting_options(vm, "all", "fixed");
        conflicting_options(vm, "all", "major-body");
        conflicting_options(vm, "all", "input");
        conflicting_options(vm, "all", "eph-file");
        conflicting_options(vm, "all", "orbit-file");
        conflicting_options(vm, "all", "horizons");
        conflicting_options(vm, "eph-file", "horizons");
        conflicting_options(vm, "eph-file", "observer");
        conflicting_options(vm, "eph-file", "fixed-target");
        conflicting_options(vm, "eph-file", "input");
        conflicting_options(vm, "input", "eph-file");
        conflicting_options(vm, "input", "orbit-file");
        conflicting_options(vm, "fixed-target", "horizons");
        conflicting_options(vm, "fixed-target", "parallax");
        conflicting_options(vm, "fixed-target", "use-uncertainty");
        conflicting_options(vm, "fixed-target", "observer");
        option_dependency(vm, "orbit-file", "horizons");

        return args;
    }
}