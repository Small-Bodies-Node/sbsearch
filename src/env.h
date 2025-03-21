#ifndef ENV_H_
#define ENV_H_

#include <cstdlib>
#include <memory>
#include <string>
#include <optional>

using std::optional;
using std::string;

namespace sbsearch
{
    struct Environment
    {
        optional<string> database;
        optional<string> log_file;
        optional<bool> verbose;

        Environment()
        {
            auto str = std::make_unique<string>();
            std::getenv("SBS_DATABASE");
        };
    };
}
#endif // ENV_H_