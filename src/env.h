#ifndef ENV_H_
#define ENV_H_

#include <string>
#include <optional>

using std::optional;
using std::string;

namespace sbsearch
{
    struct Environment
    {
        optional<string> database = {};
        optional<string> log_file = {};

        /**
         * @brief Environment variable manager.
         *
         * Reads parameters from the shell's environment, or from a .env file,
         * the latter taking precedence.
         *
         * Supported variables:
         *   - SBS_DATABASE: the database URL
         *   - SBS_LOG_FILE: a file name to use for logging
         *
         * .env file format
         *   - NAME=value
         *   - NAME is any alphanumeric character or an underscore.
         *   - Leading and trailing whitespace, and whitespace surrounding "="
         *     is ignored.
         *
         */
        Environment();
    };

    static Environment ENV;
}
#endif // ENV_H_