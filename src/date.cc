#include "config.h"

#include <iostream>
#include <string>
#include <boost/program_options.hpp>

#include "date.h"
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
            int y, m, d;
            y = (int)std::stoul(s.substr(0, 4).c_str());
            m = (int)std::stoul(s.substr(5, 2).c_str());
            d = (int)std::stoul(s.substr(8, 2).c_str());

            double djm0, djm;
            int status = iauCal2jd(y, m, d, &djm0, &mjd_);
            if (status == -1)
                throw std::range_error("Invalid year.");
            else if (status == -2)
                throw std::range_error("Invalid month.");
            else if (status == -3)
                throw std::range_error("Invalid day.");
            else if (status != 0)
                throw std::runtime_error("Unexpected status from calendar conversion.");
        }
    }

    Date::Date(const double &mjd)
    {
        mjd_ = mjd;
        int y, m, d;
        double fd;
        int status = iauJd2cal(2400000.5, mjd, &y, &m, &d, &fd);
        if (status)
            throw std::range_error("Invalid modified Julian date.");
        char buf[20];
        sprintf(buf, "%d-%02d-%02d", y, m, d);
        iso_ = string(buf);
    };

    const string Date::iso() const
    {
        return iso_;
    }

    const double Date::mjd() const
    {
        return mjd_;
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
}
