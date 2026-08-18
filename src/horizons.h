#ifndef HORIZONS_H_
#define HORIZONS_H_

#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <boost/program_options.hpp>
#include <curl/curl.h>

#include "date.h"
#include "ephemeris/ephemeris.h"
#include "moving_target.h"
#include "orbital_elements.h"

using sbsearch::ephemeris::Ephemeris;
using std::string;
using std::string_view;

namespace sbsearch
{
    // Horizons ephemeris generation.
    class Horizons
    {
    public:
        Horizons(const MovingTarget target,
                 string_view center,
                 const Date start_date,
                 const Date stop_date,
                 string_view step_size,
                 const bool cache = true);

        // Get/set properties.
        inline const MovingTarget target() { return target_; }
        inline void target(MovingTarget new_target) { target_ = new_target; }

        inline string_view center() { return center_; }
        inline void center(string new_center) { center_ = new_center; }

        inline const Date start_date() { return start_date_; }
        inline void start_date(Date new_start_date) { start_date_ = new_start_date; }

        inline const Date stop_date() { return stop_date_; }
        inline void stop_date(Date new_stop_date) { stop_date_ = new_stop_date; }

        inline string_view step_size() { return step_size_; }
        inline void step_size(string new_time_step) { step_size_ = new_time_step; }

        inline const bool cache() { return cache_; }
        inline void cache(bool new_cache) { cache_ = new_cache; }

        // The formatted Horizons ephemeris command.
        string_view command();

        // The formatted Horizons query string.
        string_view parameters();

        // The last query result as a string.
        inline string_view table() { return table_; }

        // The ephemeris data from the last query result.
        inline const Ephemeris::Data data() { return data_; }

        // Format a Horizons COMMAND.
        //
        //   target: MovingTarget
        //   mjd: set the comet closest apparition parameter to this date
        static string format_command(const MovingTarget &target, const double mjd = 0);

        // format COMMAND for specific target types
        static string major_body_command(const MovingTarget &target);
        static string small_body_command(const MovingTarget &target, const double mjd = 0);
        static string asteroid_command(const MovingTarget &target);
        static string comet_command(const MovingTarget &target, const double mjd = 0);
        static string orbit_command(const MovingTarget &target);

        // Format and store the Horizons COMMAND.
        void format_command();

        // Format a query string.
        static string format_query(string_view command,
                                   string_view center,
                                   const Date &start_date,
                                   const Date &stop_date,
                                   string_view step_size);

        // Format and store the Horizons query string.
        void format_query();

        // Get a query, possibly cached, from Horizons as a string.
        static string query(string_view parameters, const bool cache = true);

        // Get and store the Horizons query, possibly using the cache.
        void query();

        // Parse a Horizons query result (e.g., from a cached file) into an
        // ephemeris object.
        static Ephemeris::Data parse(string_view table);

        // Parse the stored Horizons query result into an ephemeris object.
        void parse();

        // Run the Horizons query and return ephemeris data.
        Ephemeris::Data get_ephemeris_data();

        // wrap s with COMMAND='s'
        static string command(string_view s);

    private:
        bool cache_;
        string center_, step_size_, command_, parameters_, table_;
        Date start_date_, stop_date_;
        MovingTarget target_;
        Ephemeris::Data data_;
    };
}

#endif // HORIZONS_H_