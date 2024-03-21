#ifndef HORIZONS_H_
#define HORIZONS_H_

#include <iostream>
#include <stdexcept>
#include <string>
#include <boost/program_options.hpp>
#include <curl/curl.h>

#include "ephemeris.h"
#include "moving_target.h"

namespace sbsearch
{
    namespace horizons
    {
        // Format a target for the Horizons COMMAND.
        const string format_command(const string target,
                                    const double mjd = 0,
                                    const bool small_body = true);

        // Format Horizons query string.
        const string format_query(const string command,
                                  const string center,
                                  const Date start_date,
                                  const Date stop_date,
                                  const string time_step);

        // Get a query from Horizons as a string.
        const string query(const string parameters, const bool cache = true);

        // Parse a Horizons query into an ephemeris object.
        Ephemeris::Data parse(const string &table);
    }
}

#endif // HORIZONS_H_