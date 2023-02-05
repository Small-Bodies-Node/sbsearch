#ifndef INDEXER_H_
#define INDEXER_H_

#include <string>
#include <vector>
#include <s2/s2region.h>
#include <s2/s2region_term_indexer.h>

#include "ephemeris.h"
#include "observation.h"

using std::string;
using std::vector;

namespace sbsearch
{
    class Indexer
    {
    public:
        // indexer options
        class Options
        {
        public:
            void max_spatial_cells(const int n);
            int max_spatial_cells() const;

            int max_spatial_level() const;
            int min_spatial_level() const;

            double max_spatial_resolution() const;
            double min_spatial_resolution() const;
            int temporal_resolution() const;

            void max_spatial_resolution(const double radians);
            void min_spatial_resolution(const double radians);

            // temporal resolution
            void temporal_resolution(const int inverse_days);

        private:
            int max_spatial_cells_ = 8;
            int min_spatial_level_ = 4;
            int max_spatial_level_ = 12;
            int time_terms_per_day_ = 100;
        };

        // Constructs an Indexer with the given Options.
        Indexer(const Options &options);

        // with the default Options.
        Indexer() : Indexer(Options()){};

        // access options
        const Options &options();

        // spatial-only index for a point
        vector<string> index_terms(const S2Point &point);
        vector<string> query_terms(const S2Point &point);

        // spatial-only index for a region
        vector<string> index_terms(const S2Region &region);
        vector<string> query_terms(const S2Region &region);

        // index the region over a time period
        vector<string> index_terms(const S2Region &region, double mjd_start, double mjd_stop);
        vector<string> query_terms(const S2Region &region, double mjd_start, double mjd_stop);

        // higher-level object indexing
        vector<string> index_terms(const Observation &observation);
        vector<string> query_terms(const Observation &observation);
        vector<string> index_terms(const Ephemeris &eph);
        vector<string> query_terms(const Ephemeris &eph);

    private:
        Options options_;
        S2RegionTermIndexer indexer_;

        enum TermStyle
        {
            index,
            query
        };

        vector<string> temporal_terms(const double mjd_start, const double mjd_stop);
        vector<string> generate_terms(const TermStyle style, const S2Point &point);
        vector<string> generate_terms(const TermStyle style, const S2Region &region);
        vector<string> generate_terms(const TermStyle style, const S2Region &region, double mjd_start, double mjd_stop);
    };
}

#endif // INDEXER_H_
