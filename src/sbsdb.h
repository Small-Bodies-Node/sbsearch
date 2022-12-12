#ifndef SBSDB_H_
#define SBSDB_H_

#include "util.h"
#include "ephemeris.h"
#include "observation.h"

#include <s2/s2point.h>
#include <s2/s2cap.h>
#include <s2/s2metrics.h>
#include <s2/s2polyline.h>
#include <s2/s2polygon.h>
#include <s2/s2region_term_indexer.h>

namespace sbsearch
{
    struct Found
    {
        Observation obs;
        Ephemeris eph;
    };

    class SBSearchDatabase
    {
    public:
        class Options
        {
        public:
            void max_spatial_cells(int n);
            int max_spatial_cells() const;

            int max_spatial_level() const;
            int min_spatial_level() const;

            double max_spatial_resolution() const;
            double min_spatial_resolution() const;
            double temporal_resolution() const;

            void max_spatial_resolution(double arcmin);
            void min_spatial_resolution(double arcmin);

            // temporal resolution (nearest integer fraction between 0.01 and 1)
            void temporal_resolution(double days);

        private:
            int max_spatial_cells_ = 8;
            int min_spatial_level_ = 4;
            int max_spatial_level_ = 12;
            int time_terms_per_day_ = 100;
        };

        // Constructs a database with the given options.
        SBSearchDatabase(const Options &options);

        // Default constructor.  Options can be set using mutable_options().
        SBSearchDatabase() : SBSearchDatabase(Options()){};

        // Returns the current options.
        const Options &options() const { return options_; }

        // close database connection
        virtual void close() = 0;

        // initialize database, or add any missing tables, indices, etc.
        virtual void setup_tables() = 0;

        void drop_time_indices();

        // void add_moving_target(const MovingTarget target);
        // MovingTarget get_moving_target(const int64 object_id);
        // MovingTarget get_moving_target(const char* name);

        // add an ephemeris to the database
        // void add_ephemeris(const Ephemeris ephemeris);
        // get an ephemeris from the database
        // Ephemeris get_ephemeris(const MovingTarget target);

        // add an observation to the database
        // - if the observation ID is not set, it will be updated
        // - if the index terms are not defined, they will be generated
        virtual void add_observation(Observation observation) = 0;

        // add a set of observations to the database, see add_observation for details
        void add_observations(vector<Observation> &observations);

        // get an observation from the database
        virtual Observation get_observation(const int64 observation_id) = 0;

        // get a set of observations from the database by observation_id, from first up to last.
        template <typename ForwardIterator>
        vector<Observation> get_observations(const ForwardIterator &first, const ForwardIterator &last);

        // approximate search functions
        virtual vector<Observation> fuzzy_search(vector<string> terms) = 0;
        virtual vector<Observation> fuzzy_search(Ephemeris eph) = 0;

        // precise search functions
        vector<Found> find_observations(S2Point point);
        vector<Found> find_observations(S2Cap cap);
        vector<Found> find_observations(S2Polyline polyline);
        vector<Found> find_observations(S2Polygon polygon);
        vector<Found> find_observations(Ephemeris ephemeris);

    protected:
        S2RegionTermIndexer indexer;

    private:
        Options options_;
        virtual void execute_sql(const char *statement) = 0;
    };

    template <typename ForwardIterator>
    vector<Observation> SBSearchDatabase::get_observations(const ForwardIterator &first, const ForwardIterator &last)
    {
        vector<Observation> observations;
        ForwardIterator observation_id = first;
        while (observation_id != last)
        {
            observations.push_back(get_observation(*observation_id));
            observation_id++;
        }
        return observations;
    }
}
#endif // SBSDB_H_
