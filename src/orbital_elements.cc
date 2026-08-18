#include <algorithm>
#include <iostream>
#include <regex>
#include <string>

#include "exceptions.h"
#include "orbital_elements.h"

using std::string;

namespace sbsearch
{
    bool operator==(const OrbitalElements &a, const OrbitalElements &b)
    {
        return (a.epoch == b.epoch) &&
               (a.ec == b.ec) &&
               (a.qr == b.qr) &&
               (a.Tp == b.Tp) &&
               (a.om == b.om) &&
               (a.w == b.w) &&
               (a.in == b.in) &&
               (a.ma == b.ma) &&
               (a.a == b.a) &&
               (a.n == b.n);
    }

    std::istream &operator>>(std::istream &is, OrbitalElements &orbit)
    {
        std::regex pattern{R"(^(\w+)=(.+)$)"};
        for (string line; getline(is, line);)
        {
            std::smatch matches;
            if (std::regex_search(line, matches, pattern))
            {
                string key = matches[1];
                std::transform(key.begin(), key.end(), key.begin(),
                               [](unsigned char c)
                               { return std::tolower(c); });

                long double value = std::stold(matches[2]);

                if (key == "epoch")
                    orbit.epoch = value;
                else if (key == "ec")
                    orbit.ec = value;
                else if (key == "qr")
                    orbit.qr = value;
                else if (key == "tp")
                    orbit.Tp = value;
                else if (key == "om")
                    orbit.om = value;
                else if (key == "w")
                    orbit.w = value;
                else if (key == "in")
                    orbit.in = value;
                else if (key == "ma")
                    orbit.ma = value;
                else if (key == "a")
                    orbit.a = value;
                else if (key == "n")
                    orbit.n = value;
                else
                    throw OrbitError("Invalid orbit parameter input \"" + line + '"');
            }
        }
        return is;
    }
}