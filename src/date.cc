#include "config.h"

#include <ctime>
#include <iostream>
#include <optional>
#include <string>
#include <boost/program_options.hpp>

#include "date.h"
#include "util/string.h"
#include "sofa/sofa.h"

using std::string;
using std::to_string;

namespace sbsearch
{
    Date::Date(const string &s)
    {
        bool iso_like = false;

        // At least 10 characters and two hyphens?  Probably a ISO formatted date.
        if ((s.size() >= 10) and std::count(s.begin(), s.end(), '-'))
            iso_like = true;

        if (!iso_like)
        {
            // Maybe it is a number for MJD?
            double mjd;
            try
            {
                mjd = std::stod(s);
            }
            catch (std::invalid_argument)
            {
                throw std::invalid_argument("String does not look like YYYY-MM-DD format and failed parsing as MJD.");
            }

            Date d(mjd);
            mjd_ = mjd;
            iso_ = d.iso();
        }
        else
        {
            iso_ = s;

            // Parse and convert the date to MJD (double).
            auto parts = util::split(s, ' ');
            int y, m, d;
            y = (int)std::stoul(parts[0].substr(0, 4).c_str());
            m = (int)std::stoul(parts[0].substr(5, 2).c_str());
            d = (int)std::stoul(parts[0].substr(8, 2).c_str());

            double fracday = 0;
            if (parts.size() > 1)
            {
                auto hms = util::split(parts[1], ':');
                if (hms.size() > 0)
                    fracday = std::stod(hms[0]) / 24.;
                if (hms.size() > 1)
                    fracday += std::stod(hms[1]) / 1440.;
                if (hms.size() > 2)
                    fracday += std::stod(hms[2]) / 86400.;
            }

            double djm0, djm, mjd;
            int status = iauCal2jd(y, m, d, &djm0, &mjd);
            if (status == -1)
                throw std::range_error("Invalid year.");
            else if (status == -2)
                throw std::range_error("Invalid month.");
            else if (status == -3)
                throw std::range_error("Invalid day.");
            else if (status != 0)
                throw std::runtime_error("Unexpected status from calendar conversion.");
            mjd_ = mjd + fracday;
            iso_ = Date(mjd_.value()).iso();
        }
    }

    Date::Date(const double &mjd)
    {
        mjd_ = mjd;
        int y, m, d;
        double fd;
        int status = iauJd2cal(2400000.5, mjd, &y, &m, &d, &fd);
        if (status || mjd > 2400000)
            throw std::range_error("Invalid modified Julian date.");

        char sign;
        int ihmsf[4];
        iauD2tf(0, fd, &sign, ihmsf);

        // did it round up to 24 hrs?
        if (ihmsf[0] == 24)
        {
            d += 1;
            ihmsf[0] = 0;
        }

        char buf[20];
        sprintf(buf, "%d-%02d-%02d %02d:%02d:%02d", y, m, d, ihmsf[0], ihmsf[1], ihmsf[2]);
        iso_ = string(buf);
    };

    Date Date::now()
    {
        char now[32];
        std::time_t time_now = std::time(nullptr);
        std::strftime(now, 32, "%F %T", std::gmtime(&time_now));
        return Date(string(now));
    }

    const string Date::iso() const
    {
        return iso_;
    }

    const double Date::mjd() const
    {
        return mjd_.value_or(-1);
    }

    Date operator+(const Date &a, const double days)
    {
        return Date(a.mjd() + days);
    }

    Date operator-(const Date &a, const double days)
    {
        return Date(a.mjd() - days);
    }

    bool operator==(const Date &a, const Date &b)
    {
        return a.mjd() == b.mjd();
    }

    bool operator>=(const Date &a, const Date &b)
    {
        return a.mjd() >= b.mjd();
    }

    bool operator<=(const Date &a, const Date &b)
    {
        return a.mjd() <= b.mjd();
    }

    bool operator>(const Date &a, const Date &b)
    {
        return a.mjd() > b.mjd();
    }

    bool operator<(const Date &a, const Date &b)
    {
        return a.mjd() < b.mjd();
    }

    std::ostream &operator<<(std::ostream &os, const Date &date)
    {
        os << date.iso();
        return os;
    }

    std::istream &operator>>(std::istream &is, Date &date)
    {
        string s;
        is >> s;
        date = Date(s);
        return is;
    }

    std::istream &operator>>(std::istream &in, DateFormat &date_format)
    {
        std::string token;
        in >> token;
        std::transform(token.begin(), token.end(), token.begin(),
                       [](unsigned char c)
                       { return std::tolower(c); });
        if (token == "mjd")
            date_format = DateFormat::MJD;
        else if (token == "calendar")
            date_format = DateFormat::CALENDAR;
        else
            in.setstate(std::ios_base::failbit);
        return in;
    }

}
