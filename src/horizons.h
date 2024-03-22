#ifndef HORIZONS_H_
#define HORIZONS_H_

#include <iostream>
#include <stdexcept>
#include <string>
#include <boost/program_options.hpp>
#include <curl/curl.h>

#include "date.h"
#include "ephemeris.h"
#include "moving_target.h"

namespace sbsearch
{
    class Horizons
    {
    public:
        // Initialize the command and query parameters, and the cache flag.
        Horizons(const MovingTarget target,
                 const string center,
                 const Date start_date,
                 const Date stop_date,
                 const string time_step,
                 const bool cache = true);

        // Get/set properties.
        inline const MovingTarget target() { return target_; }
        inline void target(MovingTarget new_target) { target_ = new_target; }

        inline const string center() { return center_; }
        inline void center(string new_center) { center_ = new_center; }

        inline const Date start_date() { return start_date_; }
        inline void start_date(Date new_start_date) { start_date_ = new_start_date; }

        inline const Date stop_date() { return stop_date_; }
        inline void stop_date(Date new_stop_date) { stop_date_ = new_stop_date; }

        inline const string time_step() { return time_step_; }
        inline void time_step(string new_time_step) { time_step_ = new_time_step; }

        inline const bool cache() { return cache_; }
        inline void cache(bool new_cache) { cache_ = new_cache; }

        // The formatted Horizons ephemeris command.
        inline const string command()
        {
            format_command();
            return command_;
        }

        // The formatted Horizons query string.
        inline const string parameters()
        {
            format_query();
            return parameters_;
        }

        // The last query result as a string.
        inline const string table() { return table_; }

        // The ephemeris data from the last query result.
        inline const Ephemeris::Data data() { return data_; }

        // Format a Horizons COMMAND.
        //
        //   designation: any Horizons resolvable string
        //   small_body: true if this is a small body object
        //   mjd: set the comet closest apparition parameter to this date
        static string format_command(const string designation,
                                     const bool small_body = true,
                                     const double mjd = 0);

        // Format and store the Horizons COMMAND.
        void format_command();

        // Format a query string.
        static string format_query(const string command,
                                   const string center,
                                   const Date start_date,
                                   const Date stop_date,
                                   const string time_step);

        // Format and store the Horizons query string.
        void format_query();

        // Get a query, possibly cached, from Horizons as a string.
        static string query(const string parameters, const bool cache = true);

        // Get and store the Horizons query, possibly using the cache.
        void query();

        // Parse a Horizons query result (e.g., from a cached file) into an
        // ephemeris object.
        static Ephemeris::Data parse(const string &table);

        // Parse the stored Horizons query result into an ephemeris object.
        void parse();

        // Run the Horizons query and return ephemeris data.
        Ephemeris::Data get_ephemeris();

    private:
        bool cache_;
        string center_, time_step_, command_, parameters_, table_;
        Date start_date_, stop_date_;
        MovingTarget target_;
        Ephemeris::Data data_;
    };
}

#endif // HORIZONS_H_