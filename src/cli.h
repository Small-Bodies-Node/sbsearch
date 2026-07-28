#ifndef CLI_H_
#define CLI_H_

#include <sstream>
#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "date.h"
#include "ephemeris.h"
#include "env.h"
#include "sbsearch.h"
#include "sofa/sofa.h"

namespace po = boost::program_options;
using std::string;
using std::vector;

namespace sbsearch
{
    namespace cli
    {
        // Auxiliary functions for checking input for validity. From or based on
        // libboost examples.

        // Check that 'opt1' and 'opt2' are not specified at the same time.
        void conflicting_options(const po::variables_map &vm,
                                 const char *opt1, const char *opt2);

        // Check that 'option' is not specified for 'action'.
        void action_conflicting_option(const po::variables_map &vm,
                                       const char *action, const char *opt2);

        // Check that if 'for_what' is specified, then 'required_option' is
        // specified too.
        void option_dependency(const po::variables_map &vm,
                               const char *for_what, const char *required_option);

        // Check that 'required_option' is specified for the given action.
        void action_dependency(const po::variables_map &vm,
                               const char *action, const char *required_option);

        bool confirm(std::string_view prompt);

        void message(std::string_view str);

        vector<std::pair<Date, Date>> date_ranges(const Date &start, const Date &stop, const double chunk);

        struct CommonArguments
        {
            string database = ENV.database.value_or("");
            string log_file = ENV.log_file.value_or("/dev/null");
            bool verbose = false;

            LogLevel log_level()
            {
                return verbose ? DEBUG : INFO;
            }
        };

        // Get the common options description.
        po::options_description get_common_options(CommonArguments *args);

        enum OutputFormat
        {
            CSV,   // .csv
            JSON,  // .json
            TABLE, // .anything-else
            AUTO   // determine by suffix
        };

        // Determine output format from a file name.  JSON if ending with .json,
        // TABLE otherwise.
        OutputFormat get_output_format(const string filename);

        // Get output format from a string (for command line argument parsing with BOOST)
        std::istream &operator>>(std::istream &in, sbsearch::cli::OutputFormat &format);
    }

    template <typename T, typename Values>
    void validate(boost::any &v, Values const &values, std::optional<T> *, int)
    {
        // H/T https://stackoverflow.com/questions/66539770/using-boostprogram-options-with-stdoptional
        po::validators::check_first_occurrence(v);
        v = std::make_optional(
            boost::lexical_cast<T>(
                po::validators::get_single_string(values)));
    }
}

#endif // CLI_H_