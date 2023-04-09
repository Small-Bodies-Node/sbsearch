#ifndef SBS_CLI_H_
#define SBS_CLI_H_

#include <sstream>
#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "ephemeris.h"
#include "sofa/sofa.h"

using std::string;
using std::vector;

namespace sbsearch
{
    namespace cli
    {
        // Auxiliary functions for checking input for validity. From or based on
        // libboost examples.

        // Check that 'opt1' and 'opt2' are not specified at the same time.
        void conflicting_options(const boost::program_options::variables_map &vm,
                                 const char *opt1, const char *opt2);

        // Check that 'option' is not specified for 'action'.
        void action_conflicting_option(const boost::program_options::variables_map &vm,
                                       const char *action, const char *opt2);

        // Check that if 'for_what' is specified, then 'required_option' is
        // specified too.
        void option_dependency(const boost::program_options::variables_map &vm,
                               const char *for_what, const char *required_option);

        // Check that 'required_option' is specified for the given action.
        void action_dependency(const boost::program_options::variables_map &vm,
                               const char *action, const char *required_option);

        struct Date
        {
            string ymd = "";
            double mjd = UNDEF_TIME;
        };

        // Overload the boost::json 'validate' function for dates. It makes sure
        // that the value is of form YYYY-MM-DD.
        void validate(boost::any &v,
                      const vector<string> &values,
                      Date *, int);

        bool confirm(const string prompt);

        // Read file contents into a string.
        const string read_file(const string &file);

        // Write HTTP data to a string.
        size_t save_http_data(void *buffer, size_t size, size_t nmemb, void *data);

        // Format the given target name as a Horizons query COMMAND string. If
        // the target appears to be a comet, then fragment searching is disabled
        // (NOFRAG) and the closest apparition is requested (CAP).  Possible
        // temporary asteroidal designations (e.g., 2000 XY) are prefixed with
        // "DES="".
        const string format_horizons_command(const string target);

        // Get an ephemeris table  Horizons's API.
        const string from_horizons(const string target,
                                   const string center,
                                   const Date start_date,
                                   const Date stop_date,
                                   const string time_step,
                                   const bool verbose);

        // Parse a Horizons ephemeris table.
        Ephemeris::Data parse_horizons(const string &table);
    }
}
#endif // SBS_CLI_H_