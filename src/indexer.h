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
            // The maximum number of cells to generate for indexing.
            void max_spatial_index_cells(const int n);
            int max_spatial_index_cells() const;

            // The maximum number of cells to generate for a query.
            void max_spatial_query_cells(const int n);
            int max_spatial_query_cells() const;

            // The maximum spatial level to consider.
            int max_spatial_level() const;
            void max_spatial_level(const int level);

            // The minimum spatial level to consider.
            int min_spatial_level() const;
            void min_spatial_level(const int level);

            // The maximum spatial scale to consider.
            double max_spatial_resolution() const;
            void max_spatial_resolution(const double radians);

            // The minimum spatial scale to consider.
            double min_spatial_resolution() const;
            void min_spatial_resolution(const double radians);

            // The temporal resolution.
            int temporal_resolution() const;
            void temporal_resolution(const int inverse_days);

            bool operator==(const Options &other) const;
            bool operator!=(const Options &other) const;

        private:
            int max_spatial_index_cells_ = 30;
            int max_spatial_query_cells_ = 8;
            int min_spatial_level_ = 4;
            int max_spatial_level_ = 12;
            int time_terms_per_day_ = 1;
        };

        // For mutable options, only max_spatial_query_cells is settable.
        class MutableOptions : public Options
        {
        public:
            void max_spatial_index_cells(const int n) = delete;
            void max_spatial_level(const int level) = delete;
            void min_spatial_level(const int level) = delete;
            void max_spatial_resolution(const double radians) = delete;
            void min_spatial_resolution(const double radians) = delete;
            void temporal_resolution(const int inverse_days) = delete;
        };

        // Constructs an Indexer with the given Options.
        Indexer(const Options &options);

        // with the default Options.
        Indexer() : Indexer(Options()){};

        // access options as a constant
        const Options &options();

        // mutable options
        MutableOptions &mutable_options();

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
