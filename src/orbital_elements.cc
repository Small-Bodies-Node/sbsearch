#include <algorithm>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <string_view>
#include <boost/json.hpp>

#include "exceptions.h"
#include "json.h"
#include "orbital_elements.h"
#include "util/string.h"

#define FLOAT_PATTERN R"((?:\d*\.\d\d*|\d+)(?:[eE][+-]?\d+)?)"

using sbsearch::json::get_string_as_long_double;
using sbsearch::json::get_string_as_optional_long_double;
using std::string;
using std::string_view;

namespace sbsearch
{
    bool OrbitalElements::is_complete() const
    {
        // Need one of {TP, QR}, {MA, A} or {MA, N}
        if (Tp.has_value() && qr.has_value())
            return true;
        if (ma.has_value() && a.has_value())
            return true;
        if (ma.has_value() && n.has_value())
            return true;
        return false;
    }

    OrbitalElements OrbitalElements::from_json(boost::json::object &obj)
    {
        sbsearch::OrbitalElements orbit{
            get_string_as_long_double(obj, {"epoch", "EPOCH"}),
            get_string_as_long_double(obj, {"ec", "e", "EC"}),
            get_string_as_optional_long_double(obj, {"qr", "q", "QR"}),
            get_string_as_optional_long_double(obj, {"Tp", "TP"}),
            get_string_as_long_double(obj, {"om", "node", "OM"}),
            get_string_as_long_double(obj, {"w", "peri", "W"}),
            get_string_as_long_double(obj, {"in", "i", "IN"}),
            get_string_as_optional_long_double(obj, {"ma", "m", "MA", "M"}),
            get_string_as_optional_long_double(obj, {"a", "A"}),
            get_string_as_optional_long_double(obj, {"n", "N"}),
        };

        return orbit;
    }

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

    bool operator!=(const OrbitalElements &a, const OrbitalElements &b)
    {
        return !(a == b);
    }

    std::ostream &operator<<(std::ostream &os, const OrbitalElements &orbit)
    {
        auto value_or_null = [&os](optional<long double> value)
        {
            if (value)
                os << *value;
            else
                os << "null";
        };

        os << std::setprecision(16)
           << "EPOCH=" << orbit.epoch
           << "\n"
           << "EC=" << orbit.ec
           << " QR=";
        value_or_null(orbit.qr);
        os << " TP=";
        value_or_null(orbit.Tp);
        os << "\n"
           << "OM=" << orbit.om
           << " W=" << orbit.w
           << " IN=" << orbit.in
           << "\n"
           << "MA=";
        value_or_null(orbit.ma);
        os << " A=";
        value_or_null(orbit.a);
        os << " N=";
        value_or_null(orbit.n);
        return os;
    }

    std::istream &operator>>(std::istream &is, OrbitalElements &orbit)
    {
        std::regex pattern{R"/((\w+)\s*=\s*()/" FLOAT_PATTERN "|null)"};
        for (string line; getline(is, line);)
        {
            if (line.size() == 0 || line[0] == '#')
                continue;

            auto start = std::sregex_iterator(line.begin(), line.end(), pattern);
            auto end = std::sregex_iterator();
            for (std::sregex_iterator i = start; i != end; i++)
            {
                std::smatch matches = *i;
                string key = matches[1].str();
                std::transform(key.begin(), key.end(), key.begin(),
                               [](unsigned char c)
                               { return std::tolower(c); });

                optional<long double> value;
                if (matches[2].str() != "null")
                    value = std::stold(matches[2].str());

                if (key == "epoch")
                    orbit.epoch = *value;
                else if ((key == "ec") || (key == "e"))
                    orbit.ec = *value;
                else if ((key == "qr") || (key == "q"))
                    orbit.qr = value;
                else if (key == "tp")
                    orbit.Tp = value;
                else if ((key == "om") || (key == "node"))
                    orbit.om = *value;
                else if ((key == "w") || (key == "peri"))
                    orbit.w = *value;
                else if ((key == "in") || (key == "i"))
                    orbit.in = *value;
                else if ((key == "ma") || (key == "m"))
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