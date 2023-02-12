#include "indexer.h"

#include <string>
#include <vector>
#include <s2/s2metrics.h>
#include <s2/s2region.h>

#include "ephemeris.h"
#include "observation.h"
#include "util.h"

// for testing
#include "sbsearch_testing.h"

using sbsearch::Indexer;
using std::string;
using std::vector;

void Indexer::Options::max_spatial_cells(const int n) { max_spatial_cells_ = n; };
int Indexer::Options::max_spatial_cells() const { return max_spatial_cells_; };
int Indexer::Options::max_spatial_level() const { return max_spatial_level_; };
int Indexer::Options::min_spatial_level() const { return min_spatial_level_; };

void Indexer::Options::max_spatial_resolution(const double radians)
{
    min_spatial_level_ = S2::kAvgEdge.GetClosestLevel(radians);
};

double Indexer::Options::max_spatial_resolution() const
{
    return S2::kAvgEdge.GetValue(min_spatial_level_);
};

void Indexer::Options::min_spatial_resolution(const double radians)
{
    max_spatial_level_ = S2::kAvgEdge.GetClosestLevel(radians);
};

double Indexer::Options::min_spatial_resolution() const
{
    return S2::kAvgEdge.GetValue(max_spatial_level_);
};

void Indexer::Options::temporal_resolution(const int inverse_days)
{
    time_terms_per_day_ = inverse_days;
};

int Indexer::Options::temporal_resolution() const
{
    return time_terms_per_day_;
};

Indexer::Indexer(const Options &options)
{
    options_ = options;
    S2RegionTermIndexer::Options s2options;
    s2options.set_max_cells(options.max_spatial_cells());
    s2options.set_min_level(options.min_spatial_level());
    s2options.set_max_level(options.max_spatial_level());
    indexer_ = S2RegionTermIndexer(s2options);
}

const Indexer::Options &Indexer::options()
{
    return options_;
}

vector<string> Indexer::index_terms(const S2Point &point)
{
    return generate_terms(index, point);
}

vector<string> Indexer::query_terms(const S2Point &point)
{
    return generate_terms(query, point);
}

vector<string> Indexer::index_terms(const S2Region &region)
{
    return generate_terms(index, region);
}

vector<string> Indexer::query_terms(const S2Region &region)
{
    return generate_terms(query, region);
}

vector<string> Indexer::index_terms(const S2Region &region, double mjd_start, double mjd_stop)
{
    return generate_terms(index, region, mjd_start, mjd_stop);
}

vector<string> Indexer::query_terms(const S2Region &region, double mjd_start, double mjd_stop)
{
    return generate_terms(query, region, mjd_start, mjd_stop);
}

vector<string> Indexer::index_terms(const Observation &observation)
{
    return generate_terms(index, observation.as_polygon(), observation.mjd_start(), observation.mjd_stop());
}

vector<string> Indexer::query_terms(const Observation &observation)
{
    return generate_terms(query, observation.as_polygon(), observation.mjd_start(), observation.mjd_stop());
}

vector<string> Indexer::index_terms(const Ephemeris &eph)
{
    std::set<string> all_terms;
    vector<string> segment_terms;
    for (auto segment : eph.segments())
    {
        segment_terms = generate_terms(index, segment.as_polygon(), segment.mjd(0), segment.mjd(1));
        all_terms.insert(segment_terms.begin(), segment_terms.end());
    }
    return vector<string>(all_terms.begin(), all_terms.end());
}

vector<string> Indexer::query_terms(const Ephemeris &eph)
{
    vector<string> all_terms, segment_terms;
    for (auto segment : eph.segments())
    {
        S2Polygon polygon = segment.as_polygon();
        segment_terms = generate_terms(query, polygon, segment.mjd(0), segment.mjd(1));
        // segment_terms = generate_terms(query, polygon);
        all_terms.insert(all_terms.end(), segment_terms.begin(), segment_terms.end());

        vector<S2LatLng> coords(polygon.loop(0)->vertices_span().size());
        std::transform(polygon.loop(0)->vertices_span().begin(), polygon.loop(0)->vertices_span().end(), coords.begin(),
                       [](const S2Point &p)
                       { return S2LatLng(p); });
    }
    return all_terms;
}

vector<string> Indexer::temporal_terms(const double mjd_start, const double mjd_stop)
{
    vector<string> terms;
    unsigned int left_term, right_term;
    left_term = (unsigned int)floor(mjd_start * options_.temporal_resolution());
    right_term = (unsigned int)ceil(mjd_stop * options_.temporal_resolution());

    for (unsigned int i = left_term; i < right_term; i++)
        terms.push_back(std::to_string(i));

    return terms;
}

vector<string> Indexer::generate_terms(const TermStyle style, const S2Point &point)
{
    // spatial terms
    return (style == index)
               ? indexer_.GetIndexTerms(point, "")
               : indexer_.GetQueryTerms(point, "");
}

vector<string> Indexer::generate_terms(const TermStyle style, const S2Region &region)
{
    // spatial terms
    return (style == index)
               ? indexer_.GetIndexTerms(region, "")
               : indexer_.GetQueryTerms(region, "");
}

vector<string> Indexer::generate_terms(const TermStyle style, const S2Region &region, double mjd_start, double mjd_stop)
{
    // spatial terms
    vector<string> s_terms = generate_terms(style, region);

    // temporal terms
    vector<string> t_terms = temporal_terms(mjd_start, mjd_stop);

    // Join query terms, each segment gets a time suffix, save to terms string
    vector<string> terms;
    for (auto t : t_terms)
        for (auto s : s_terms)
            terms.push_back(s + "-" + t);

    return terms;
}