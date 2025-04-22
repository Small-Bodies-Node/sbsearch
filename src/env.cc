#include <cstdlib>
#include <fstream>
#include <regex>
#include <string>
#include <optional>

#include "env.h"

using std::endl;
using std::optional;
using std::string;

namespace sbsearch
{
    Environment::Environment()
    {
        // First, get values from environment
        char *var = std::getenv("SBS_DATABASE");
        if (var)
            database.emplace(var);

        var = std::getenv("SBS_LOG_FILE");
        if (var)
            log_file.emplace(var);

        // Then, check in .env file
        std::ifstream input(".env", input.binary | input.in);
        if (input.is_open())
        {
            std::regex pattern("\\s*(\\w+)\\s*=\\s*(.+?)\\s*(?=$)");
            std::smatch matches;
            for (string line; std::getline(input, line);)
            {
                std::regex_search(line, matches, pattern);
                if (!matches.empty())
                {
                    if (matches[1] == "SBS_DATABASE")
                        database.emplace(matches[2]);
                    else if (matches[1] == "SBS_LOG_FILE")
                        log_file.emplace(matches[2]);
                }
            }
        }
    }
}