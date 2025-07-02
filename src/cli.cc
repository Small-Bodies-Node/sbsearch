#include <fstream>
#include <iostream>
#include <cstdio>
#include <streambuf>
#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "date.h"
#include "env.h"
#include "cli.h"

namespace po = boost::program_options;
using std::cerr;
using std::string;
using std::vector;

namespace sbsearch
{
    namespace cli
    {
        void conflicting_options(const po::variables_map &vm,
                                 const char *opt1, const char *opt2)
        {
            if (vm.count(opt1) && !vm[opt1].defaulted() && vm.count(opt2) && !vm[opt2].defaulted())
                throw std::logic_error(string("Conflicting options '") + opt1 + "' and '" + opt2 + "'.");
        }

        void action_conflicting_option(const po::variables_map &vm,
                                       const char *action, const char *option)
        {
            if (!vm.count("action"))
                return;

            if ((vm["action"].as<string>() == action) && vm.count(option) && !vm[option].defaulted())
                throw std::logic_error(string("Action '") + action + "' does not use option '" + option + "'.");
        }

        void option_dependency(const po::variables_map &vm,
                               const char *for_what, const char *required_option)
        {
            if (vm.count(for_what) && !vm[for_what].defaulted())
                if (vm.count(required_option) == 0 || vm[required_option].defaulted())
                    throw std::logic_error(string("Option '") + for_what + "' requires option '" + required_option + "'.");
        }

        void action_dependency(const po::variables_map &vm,
                               const char *action, const char *required_option)
        {
            if (!vm.count("action"))
                return;

            if ((vm["action"].as<string>() == action) && !vm.count(required_option))
                throw std::logic_error(string("Action '") + action + "' requires option '" + required_option + "'.");
        }

        bool confirm(std::string_view prompt)
        {
            string response;
            std::cout << prompt << " " << std::flush;
            std::cin >> response;
            return ((response[0] == 'y') || (response[1] == 'Y'));
        }

        void message(std::string_view str)
        {
            Logger::info() << str << std::endl;
            std::cout << str << std::endl;
        }

        std::istream &operator>>(std::istream &in, sbsearch::cli::OutputFormat &format)
        {
            std::string token;
            in >> token;
            std::transform(token.begin(), token.end(), token.begin(),
                           [](unsigned char c)
                           { return std::tolower(c); });
            if (token == "table")
                format = OutputFormat::TABLE;
            else if (token == "json")
                format = OutputFormat::JSON;
            else
                in.setstate(std::ios_base::failbit);
            return in;
        }

        po::options_description get_common_options(CommonArguments *args)
        {
            using namespace boost::program_options;

            options_description general("General options");
            general.add_options()(
                "database,D", value<string>(&args->database), "SBSearch database URI")(
                "log-file,L", value<string>(&args->log_file), "log file name")(
                "help,h", "display this help and exit")(
                "version", "output version information and exit")(
                "verbose,v", bool_switch(&args->verbose), "show debugging messages");

            return general;
        }
    }
}
