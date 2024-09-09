#ifndef DATE_H_
#define DATE_H_

#include "config.h"

#include <iostream>
#include <string>

#include "util.h"

using std::string;

namespace sbsearch
{
    class Date
    {
    public:
        Date() {};

        // Initialize from a string.  May be ISO format or MJD.
        Date(const string &s);

        // Initialize from modified Julian date.
        Date(const double &mjd);

        // Date in ISO, YYYY-MM-DD, format
        const string iso() const;

        // Date in MJD format.
        const double mjd() const;

        friend std::ostream &operator<<(std::ostream &os, const Date &date);
        friend std::istream &operator>>(std::istream &is, Date &date);

    private:
        std::string iso_ = "";
        double mjd_ = UNDEF_TIME;
    };
}

#endif // DATE_H_