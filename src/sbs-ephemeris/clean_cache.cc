#include <ctime>
#include <string>
#include <boost/filesystem.hpp>

#include "cli.h"
#include "files.h"

using namespace sbsearch;
using std::to_string;

namespace fs = boost::filesystem;

namespace sbsearch::sbs_ephemeris
{
    void clean_cache(double max_age)
    {
        fs::path path = get_cache_location();

        std::time_t now = std::time(nullptr);
        int checked = 0, removed = 0;
        for (fs::directory_entry &x : fs::directory_iterator(path))
        {
            std::time_t t = fs::last_write_time(x);
            if ((now - t) / 86400 > max_age)
            {
                fs::remove(x);
                removed++;
            }
            checked++;
        }

        cli::message::info(to_string(removed) + " of " + to_string(checked) +
                           " files older than " + to_string(max_age) + " days removed.");
    }
}