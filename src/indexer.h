#ifndef INDEXER_H_
#define INDEXER_H_

#include <string>
#include <vector>
#include <s2/s2region.h>
#include <s2/s2region_term_indexer.h>
#include <s2/s2shape_index_buffered_region.h>

#include "ephemeris/ephemeris.h"
#include "observation.h"

using sbsearch::ephemeris::Ephemeris;
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
            inline const int &max_spatial_index_cells() const { return max_spatial_index_cells_; };
            inline void max_spatial_index_cells(const int n) { max_spatial_index_cells_ = n; };

            // The maximum number of cells to generate for a query.
            inline const int &max_spatial_query_cells() const { return max_spatial_query_cells_; };
            inline void max_spatial_query_cells(const int n) { max_spatial_query_cells_ = n; };

            // The maximum spatial level to consider.
            inline const int &max_spatial_level() const { return max_spatial_level_; };
            inline void max_spatial_level(const int level) { max_spatial_level_ = level; };

            // The minimum spatial level to consider.
            inline const int &min_spatial_level() const { return min_spatial_level_; };
            inline void min_spatial_level(const int level) { min_spatial_level_ = level; };

            // The maximum spatial scale to consider.
            double max_spatial_resolution() const;
            void max_spatial_resolution(const double radians);

            // The minimum spatial scale to consider.
            double min_spatial_resolution() const;
            void min_spatial_resolution(const double radians);

            bool operator==(const Options &other) const;
            bool operator!=(const Options &other) const;

            // Some simple checks on the parameters.
            bool verify();

        private:
            int max_spatial_index_cells_ = 30;
            int max_spatial_query_cells_ = 8;
            int min_spatial_level_ = 4;
            int max_spatial_level_ = 12;
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
        };

        // Constructs an Indexer with the given Options.
        Indexer(const Options &options);

        // with the default Options.
        Indexer() : Indexer(Options()) {};

        // access options as a constant
        const Options &options();

        // mutable options
        MutableOptions &mutable_options();

        enum TermStyle
        {
            index,
            query
        };

        // spatial-only index for a point
        vector<string> terms(const TermStyle style, const S2Point &point);
        // vector<string> query_terms(const S2Point &point);

        // spatial-only index for a region
        vector<string> terms(const TermStyle style, const S2Region &region);

        // higher-level object indexing
        vector<string> terms(const TermStyle style, const Observation &observation);
        // padding in arcsec
        vector<string> terms(const TermStyle style,
                             const Ephemeris &eph,
                             bool use_uncertainty = false,
                             double padding = 0);

    private:
        Options options_;
        S2RegionTermIndexer indexer_;

        vector<string> generate_terms(const TermStyle style, const S2Point &point);
        vector<string> generate_terms(const TermStyle style, const S2Region &region);
    };
}

#endif // INDEXER_H_
