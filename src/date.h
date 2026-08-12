#ifndef DATE_H_
#define DATE_H_

#include "config.h"

#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

using std::string;
using std::string_view;
using std::vector;

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
        // Date format: modified Julian date or calendar.
        enum Format
        {
            MJD,
            CALENDAR
        };

        Date() {};

        // Initialize from a string.  May be ISO format or MJD.  YYYY-MM-DD
        // required.  hh:mm:ss.sss is optional.
        Date(string_view s);

        // Initialize from modified Julian date, must be < 2,400,000
        Date(const double &mjd);

        // New date for the current UTC time.
        static Date now();

        // Date in ISO, YYYY-MM-DD hh:mm:ss, format.  This will always be a
        // formatted version of mjd().  Fractional seconds may be truncated.
        string iso() const;

        // Date in MJD format.
        double mjd() const;

        friend Date operator+(const Date &a, const double days);
        friend Date operator-(const Date &a, const double days);
        friend bool operator==(const Date &a, const Date &b);
        friend bool operator>=(const Date &a, const Date &b);
        friend bool operator<=(const Date &a, const Date &b);
        friend bool operator>(const Date &a, const Date &b);
        friend bool operator<(const Date &a, const Date &b);
        friend std::ostream &operator<<(std::ostream &os, const Date &date);
        friend std::istream &operator>>(std::istream &is, Date &date);

    private:
        string iso_ = "";
        std::optional<double> mjd_;
    };

    // Return a set of date ranges of length `chunk` in days.
    vector<std::pair<Date, Date>> date_ranges(const Date &start, const Date &stop, const double chunk);

    std::istream &operator>>(std::istream &in, Date::Format &date_format);
}

#endif // DATE_H_