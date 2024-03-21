#ifndef SBS_CLI_H_
#define SBS_CLI_H_

#include <sstream>
#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "date.h"
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

        bool confirm(const string prompt);
    }
}

#endif // SBS_CLI_H_