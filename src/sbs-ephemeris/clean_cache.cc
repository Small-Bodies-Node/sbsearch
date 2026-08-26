#include <chrono>
#include <ctime>
#include <filesystem>
#include <string>

#include "cli.h"
#include "files.h"

using namespace sbsearch;
using namespace std::chrono_literals;
using std::to_string;

namespace fs = std::filesystem;

namespace sbsearch::sbs_ephemeris
{
    void clean_cache(double max_age)
    {
        fs::path path = get_cache_location();

        fs::file_time_type now = fs::file_time_type::clock::now();
        int checked = 0, removed = 0;
        for (auto const &f : fs::directory_iterator(path))
        {
            if (!fs::is_regular_file(f))
                continue;

            if ((now - fs::last_write_time(f)) > (max_age * 24 * 1h))
            {
                fs::remove(f);
                removed++;
            }
            checked++;
        }

        cli::message::info(to_string(removed) + " of " + to_string(checked) +
                           " files older than " + to_string(max_age) + " days removed.");
    }
}