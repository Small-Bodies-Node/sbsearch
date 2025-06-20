#ifndef DATE_H_
#define DATE_H_

#include "config.h"

#include <iostream>
#include <optional>
#include <string>

using std::string;

namespace sbsearch
{
    /**
     * @brief Date conversions.
     *
     * Convert between MJD and ISO calendar format.
     *
     * Calendar format is YYYY-MM-DD hh:mm:ss.  Seconds are truncated, but MJD
     * preserves the full precision of the input.
     *
     */
    class Date
    {
    public:
        Date() {};

        // Initialize from a string.  May be ISO format or MJD.  YYYY-MM-DD
        // required.  hh:mm:ss.sss is optional.
        Date(const string &s);

        // Initialize from modified Julian date.
        Date(const double &mjd);

        // Date in ISO, YYYY-MM-DD hh:mm:ss, format.  This will always be a
        // formatted version of mjd().  Fractional seconds are truncated.
        const string iso() const;

        // Date in MJD format.
        const double mjd() const;

        friend std::ostream &operator<<(std::ostream &os, const Date &date);
        friend std::istream &operator>>(std::istream &is, Date &date);

    private:
        string iso_ = "";
        std::optional<double> mjd_;
    };
}

#endif // DATE_H_