#ifndef SBS_FOUND_H_
#define SBS_FOUND_H_

#include "config.h"

#include <iterator>
#include <optional>
#include <string>
#include <vector>

#include "ephemeris/ephemeris.h"
#include "observation.h"

using sbsearch::ephemeris::Ephemeris;
using std::optional;
using std::string;
using std::vector;

namespace sbsearch
{
    // A moving target ephemeris and the observation it was found in.
    //
    // The ephemeris may be a segment or a single point.
    struct Found
    {
        Observation observation;
        Ephemeris ephemeris;
        double mjd_added = 0;

        Found(Observation o, Ephemeris e, double added = 0) : observation(o), ephemeris(e), mjd_added(added) {};

        // Copy constructor.
        Found(const Found &other) = default;

        // Move constructor.
        Found(Found &&other) = default;

        // Copy assignment.
        Found &operator=(const Found &other) = default;

        // Move assignment.
        Found &operator=(Found &&other) = default;

        // Compares the observation and ephemeris objects, does not compare
        // mjd_added.
        const bool operator==(const Found &other) const;
        const bool operator!=(const Found &other) const;

        // A row in the found database table.
        struct DBModel
        {
            int64_t found_id;
            int64_t observation_id;
            int64_t moving_target_id;
            double mjd;
            double tmtp;
            double ra;
            double dec;
            double mu;
            double mu_theta;
            double unc_a;
            double unc_b;
            double unc_theta;
            double rh;
            double delta;
            double phase;
            double selong;
            double true_anomaly;
            double sangle;
            double vangle;
            optional<double> vmag;
            double mjd_added;
        };
    };

    struct Founds
    {
        vector<Found> data;
        Observation::Format observation_format;
        Ephemeris::Format ephemeris_format;

        // Default constructor is an empty vector.
        Founds() {};

        // Initialize with a vector of Found.
        Founds(const vector<Found> &founds);

        // Initialize with a vector of Found.
        Founds(vector<Found> &&founds);

        // Copy constructor.
        Founds(const Founds &founds) = default;

        // Move constructor.
        Founds(Founds &&founds) = default;

        // Copy assignment
        Founds &operator=(const Founds &other) = default;

        // Move assignment
        Founds &operator=(Founds &&other) = default;

        // Access element by index.
        inline const Found &operator[](int i) const { return data[i]; };

        // Append a single found object.
        void append(const Found &found);
        void append(Found &&found);

        // Append a vector of found objects.
        void append(const vector<Found> &founds);
        void append(vector<Found> &&founds);

        // Append another Founds object.
        void append(const Founds &founds);
        void append(Founds &&founds);

        // Pointer to beginning of Found vector.
        auto begin() const { return data.begin(); };

        // Pointer to end of Found vector.
        auto end() const { return data.end(); };

        // Number of found items.
        size_t size() const { return data.size(); }

        // Data sources.
        vector<string> source() const;

        // Observatories.
        vector<string> observatory() const;

        // Observational product IDs.
        vector<string> product_id() const;

        // Observational fields-of-view.
        vector<string> fov() const;

        // Observation IDs.
        vector<optional<int64_t>> observation_id() const;

        // Observation start dates.
        vector<double> mjd_start() const;

        // Observation stop dates.
        vector<double> mjd_stop() const;

        // Observation exposure times (s).
        vector<double> exposure() const;

        // Found moving target IDs.
        vector<optional<int64_t>> moving_target_id() const;

        // Found moving target designations.
        vector<string> designation() const;

        // Found moving target small-body flags.
        vector<bool> small_body() const;

        // Times at which the ephemerides are calculated.
        vector<optional<double>> mjd() const;

        // mjd() in calendar form
        vector<string> date() const;

        // Approximate time relative to perihelion (days).
        vector<optional<double>> tmtp() const;

        // Ephemeris right ascension.
        vector<optional<double>> ra() const;

        // Ephemeris declination.
        vector<optional<double>> dec() const;

        // Ephemeris uncertainty ellipse semi-major axis.
        vector<optional<double>> unc_a() const;

        // Ephemeris uncertainty ellipse semi-minor axis.
        vector<optional<double>> unc_b() const;

        // Ephemeris proper motion.
        vector<optional<double>> mu() const;

        // Ephemeris proper motion direction.
        vector<optional<double>> mu_theta() const;

        // Ephemeris uncertainty ellipse semi-major axis position angle (deg E of N).
        vector<optional<double>> unc_theta() const;

        // Heliocentric distance (au).
        vector<optional<double>> rh() const;

        // Observer-target distance (au).
        vector<optional<double>> delta() const;

        // Sun-target-observer angle (deg).
        vector<optional<double>> phase() const;

        // Solar elongation (deg).
        vector<optional<double>> selong() const;

        // Orbital true anomaly (deg).
        vector<optional<double>> true_anomaly() const;

        // Projected target-sun vector position angle (deg).
        vector<optional<double>> sangle() const;

        // Projected target velocity vector position angle (deg).
        vector<optional<double>> vangle() const;

        // Target apparent magnitude (meaning varies depending on ephemeris source).
        vector<optional<double>> vmag() const;

        // Date found row was added to the database.
        vector<double> mjd_added() const;
    };

    std::ostream &operator<<(std::ostream &os, const Found &found);
    std::ostream &operator<<(std::ostream &os, const Founds &founds);
}

#endif // SBS_FOUND_H_
