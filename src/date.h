#ifndef DATE_H_
#define DATE_H_

#include "config.h"

#include <iostream>
#include <string>
#include <boost/program_options.hpp>

#include "util.h"

using std::string;

namespace sbsearch
{
    class Date
    {
    public:
        Date(){};

        // Initialize from a string.
        Date(const string &iso);

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

// // Overload the boost 'validate' function for dates. It makes sure
// // that the value is of form YYYY-MM-DD.
// namespace boost
// {
//     namespace program_options
//     {
//         void validate(boost::any &v, const vector<string> &values, sbsearch::Date *, int);
//     }
// }

#endif // DATE_H_