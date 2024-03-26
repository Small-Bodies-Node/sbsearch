#include "indexer.h"

#include <set>
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

void Indexer::Options::max_spatial_index_cells(const int n) { max_spatial_index_cells_ = n; };
int Indexer::Options::max_spatial_index_cells() const { return max_spatial_index_cells_; };

void Indexer::Options::max_spatial_query_cells(const int n) { max_spatial_query_cells_ = n; };
int Indexer::Options::max_spatial_query_cells() const { return max_spatial_query_cells_; };

int Indexer::Options::max_spatial_level() const { return max_spatial_level_; };
int Indexer::Options::min_spatial_level() const { return min_spatial_level_; };

void Indexer::Options::max_spatial_level(const int level) { max_spatial_level_ = level; };
void Indexer::Options::min_spatial_level(const int level) { min_spatial_level_ = level; };

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

bool Indexer::Options::operator==(const Options &other) const
{
    return ((max_spatial_level() == other.max_spatial_level()) &
            (min_spatial_level() == other.min_spatial_level()) &
            (max_spatial_index_cells() == other.max_spatial_index_cells()) &
            (max_spatial_query_cells() == other.max_spatial_query_cells()) &
            (temporal_resolution() == other.temporal_resolution()));
}

bool Indexer::Options::operator!=(const Options &other) const
{
    return !((*this) == other);
}

Indexer::Indexer(const Options &options)
{
    options_ = options;
    S2RegionTermIndexer::Options s2options;
    s2options.set_min_level(options.min_spatial_level());
    s2options.set_max_level(options.max_spatial_level());
    indexer_ = S2RegionTermIndexer(s2options);
}

const Indexer::Options &Indexer::options()
{
    return options_;
}

Indexer::MutableOptions &Indexer::mutable_options()
{
    return static_cast<MutableOptions &>(options_);
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
    S2Polygon polygon;
    observation.as_polygon(polygon);
    return generate_terms(index, polygon, observation.mjd_start(), observation.mjd_stop());
}

vector<string> Indexer::query_terms(const Observation &observation)
{
    S2Polygon polygon;
    observation.as_polygon(polygon);
    return generate_terms(query, polygon, observation.mjd_start(), observation.mjd_stop());
}

vector<string> Indexer::index_terms(const Ephemeris &eph)
{
    std::set<string> all_terms;
    vector<string> segment_terms;
    S2Polygon polygon;
    for (auto segment : eph.segments())
    {
        segment.as_polygon(polygon);
        segment_terms = generate_terms(index, polygon, segment.data(0).mjd, segment.data(1).mjd);
        all_terms.insert(segment_terms.begin(), segment_terms.end());
    }
    return vector<string>(all_terms.begin(), all_terms.end());
}

vector<string> Indexer::query_terms(const Ephemeris &eph)
{
    std::set<string> all_terms;
    vector<string> segment_terms;
    S2Polygon polygon;
    for (auto segment : eph.segments())
    {
        segment.as_polygon(polygon);
        segment_terms = generate_terms(query, polygon, segment.data(0).mjd, segment.data(1).mjd);
        all_terms.insert(segment_terms.begin(), segment_terms.end());
    }
    return vector<string>(all_terms.begin(), all_terms.end());
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
    indexer_.mutable_options()->set_max_cells((style == index)
                                                  ? options_.max_spatial_index_cells()
                                                  : options_.max_spatial_query_cells());
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