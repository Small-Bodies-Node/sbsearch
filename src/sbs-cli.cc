#include <fstream>
#include <iostream>
#include <cstdio>
#include <streambuf>
#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "date.h"
#include "sbs-cli.h"

using std::cerr;
using std::string;
using std::vector;

namespace sbsearch
{
    namespace cli
    {
        void conflicting_options(const boost::program_options::variables_map &vm,
                                 const char *opt1, const char *opt2)
        {
            if (vm.count(opt1) && !vm[opt1].defaulted() && vm.count(opt2) && !vm[opt2].defaulted())
                throw std::logic_error(string("Conflicting options '") + opt1 + "' and '" + opt2 + "'.");
        }

        void action_conflicting_option(const boost::program_options::variables_map &vm,
                                       const char *action, const char *option)
        {
            if (!vm.count("action"))
                return;

            if ((vm["action"].as<string>() == action) & vm.count(option))
                throw std::logic_error(string("Action '") + action + "' does not use option '" + option + "'.");
        }

        void option_dependency(const boost::program_options::variables_map &vm,
                               const char *for_what, const char *required_option)
        {
            if (vm.count(for_what) && !vm[for_what].defaulted())
                if (vm.count(required_option) == 0 || vm[required_option].defaulted())
                    throw std::logic_error(string("Option '") + for_what + "' requires option '" + required_option + "'.");
        }

        void action_dependency(const boost::program_options::variables_map &vm,
                               const char *action, const char *required_option)
        {
            if (!vm.count("action"))
                return;

            if ((vm["action"].as<string>() == action) & !vm.count(required_option))
                throw std::logic_error(string("Action '") + action + "' requires option '" + required_option + "'.");
        }

        bool confirm(const string prompt)
        {
            string response;
            std::cout << prompt << " " << std::flush;
            std::cin >> response;
            return ((response[0] == 'y') | (response[1] == 'Y'));
        }
    }
}
