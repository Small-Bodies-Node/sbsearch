#ifndef SBS_FOUND_H_
#define SBS_FOUND_H_

#include "config.h"

#include <iterator>
#include <string>
#include <vector>

#include "ephemeris.h"
#include "observation.h"

using std::string;
using std::vector;

namespace sbsearch
{
    // A moving target ephemeris and the observation it was found in.
    //
    // The ephemeris may be a segment or a singe point.
    struct Found
    {
        Observation observation;
        Ephemeris ephemeris;

        Found(Observation o, Ephemeris e) : observation(o), ephemeris(e){};

        const bool operator==(const Found &other) const;

        const bool operator!=(const Found &other) const;
    };

    struct Founds
    {
        vector<Found> data;

        // Default constructor is an empty vector.
        Founds(){};

        // Initialize with a vector of Found.
        Founds(const vector<Found> &founds);

        // Copy constructor.
        Founds(const Founds &founds);

        // Access element by index.
        inline const Found &operator[](int i) const { return data[i]; };

        // Append a single found object.
        void append(const Found &found);

        // Append a vector of found objects.
        void append(const vector<Found> &founds);

        // Append another Founds object.
        void append(const Founds &founds);

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
        vector<int64> observation_id() const;

        // Observation start dates.
        vector<double> mjd_start() const;

        // Observation stop dates.
        vector<double> mjd_stop() const;

        // Observation exposure times (s).
        vector<double> exposure() const;

        // Found moving target IDs.
        vector<int64> moving_target_id() const;

        // Found moving target designations.
        vector<string> designation() const;

        // Found moving target small-body flags.
        vector<bool> small_body() const;

        // Times at which the ephemerides are calculated.
        vector<double> mjd() const;

        // Approximate time relative to perihelion (days).
        vector<double> tmtp() const;

        // Ephemeris right ascension.
        vector<double> ra() const;

        // Ephemeris declination.
        vector<double> dec() const;

        // Ephemeris uncertainty ellipse semi-major axis.
        vector<double> unc_a() const;

        // Ephemeris uncertainty ellipse semi-minor axis.
        vector<double> unc_b() const;

        // Ephemeris uncertainty ellipse semi-major axis position angle (deg E of N).
        vector<double> unc_theta() const;

        // Heliocentric distance (au).
        vector<double> rh() const;

        // Observer-target distance (au).
        vector<double> delta() const;

        // Sun-target-observer angle (deg).
        vector<double> phase() const;

        // Solar elongation (deg).
        vector<double> selong() const;

        // Orbital true anomaly (deg).
        vector<double> true_anomaly() const;

        // Projected target-sun vector position angle (deg).
        vector<double> sangle() const;

        // Projected target velocity vector position angle (deg).
        vector<double> vangle() const;

        // Target apparent magnitude (meaning varies depending on ephemeris source).
        vector<double> vmag() const;
    };

    // show_fov is considered true if set on any observation
    std::ostream &operator<<(std::ostream &os, const Found &found);
    std::ostream &operator<<(std::ostream &os, const Founds &founds);
}

#endif // SBS_FOUND_H_
