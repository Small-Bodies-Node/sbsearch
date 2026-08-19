#include <charconv>
#include <string>
#include <string_view>
#include <vector>
#include <s2/s2latlng.h>
#include <s2/s2loop.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>

#include "config.h"
#include "string.h"

using std::string;
using std::string_view;
using std::vector;

namespace sbsearch::util
{
    string pluralize(int count, string_view singular, string_view plural)
    {
        return string(count == 1 ? singular : plural);
    }

    vector<string> split(std::string_view str, const char delimiter)
    {
        int start = 0, end;
        vector<string> parts;
        while ((end = str.find(delimiter, start)) != string::npos)
        {
            parts.emplace_back(str.substr(start, end - start));
            start = end + 1;
        }

        if (start < str.length())
            parts.emplace_back(str.substr(start)); // remainder of string

        return parts;
    }

    string strip(std::string_view str)
    {
        string::size_type start = str.find_first_not_of(" ");
        string::size_type stop = str.find_last_not_of(" ");
        if (stop == string::npos)
            return "";
        return string(str.substr(start, stop - start + 1));
    }

    double svtod(string_view s)
    {
        double d;
        auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), d);
        if (ec != std::errc())
            throw std::invalid_argument("String does not look like a double.");
        return d;
    }
}