#ifndef SBS_CLI_H_
#define SBS_CLI_H_

#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "sofa/sofa.h"

using std::string;
using std::vector;

// Auxiliary functions for checking input for validity. From libboost examples.

// Check that 'opt1' and 'opt2' are not specified at the same time.
void conflicting_options(const boost::program_options::variables_map &vm,
                         const char *opt1, const char *opt2)
{
    if (vm.count(opt1) && !vm[opt1].defaulted() && vm.count(opt2) && !vm[opt2].defaulted())
        throw std::logic_error(string("Conflicting options '") + opt1 + "' and '" + opt2 + "'.");
}

// Function used to check that of 'for_what' is specified, then
// 'required_option' is specified too.
void option_dependency(const boost::program_options::variables_map &vm,
                       const char *for_what, const char *required_option)
{
    if (vm.count(for_what) && !vm[for_what].defaulted())
        if (vm.count(required_option) == 0 || vm[required_option].defaulted())
            throw std::logic_error(string("Option '") + for_what + "' requires option '" + required_option + "'.");
}

struct Date
{
    string ymd = "";
    double mjd = -1;
};

// Overload the 'validate' function for dates.
// It makes sure that value is of form YYYY-MM-DD.
void validate(boost::any &v,
              const vector<string> &values,
              Date *, int)
{
    using namespace boost::program_options;

    // Make sure no previous assignment was made.
    validators::check_first_occurrence(v);

    // Extract the first string from 'values'. If there is more than
    // one string, it's an error, and exception will be thrown.
    const string &s = validators::get_single_string(values);
    if (s.size() != 10)
        throw validation_error(validation_error::invalid_option_value);

    // Pattern match and convert the date to MJD (double).
    int y, m, d;
    y = (int)std::stoul(s.substr(0, 4).c_str());
    m = (int)std::stoul(s.substr(5, 2).c_str());
    d = (int)std::stoul(s.substr(8, 2).c_str());

    double djm0, djm;
    int status = iauCal2jd(y, m, d, &djm0, &djm);
    if (status)
        throw validation_error(validation_error::invalid_option_value);

    v = boost::any(Date{s, djm});
}

#endif // SBS_CLI_H_