#include "indexer.h"

#include <set>
#include <string>
#include <vector>
#include <s2/s2metrics.h>
#include <s2/s2region.h>
#include <s2/s2lax_polyline_shape.h>
#include <s2/s2point_vector_shape.h>

#include "ephemeris.h"
#include "observation.h"
#include "util/polygon.h"

// for testing
#include "sbsearch_testing.h"

using sbsearch::Indexer;
using std::string;
using std::vector;

void Indexer::Options::max_spatial_index_cells(const int n) { max_spatial_index_cells_ = n; };
const int &Indexer::Options::max_spatial_index_cells() const { return max_spatial_index_cells_; };

void Indexer::Options::max_spatial_query_cells(const int n) { max_spatial_query_cells_ = n; };
const int &Indexer::Options::max_spatial_query_cells() const { return max_spatial_query_cells_; };

const int &Indexer::Options::max_spatial_level() const { return max_spatial_level_; };
void Indexer::Options::max_spatial_level(const int level) { max_spatial_level_ = level; };

const int &Indexer::Options::min_spatial_level() const { return min_spatial_level_; };
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

bool Indexer::Options::operator==(const Options &other) const
{
    return std::tie(max_spatial_level_,
                    min_spatial_level_,
                    max_spatial_index_cells_,
                    max_spatial_query_cells_) ==
           std::tie(other.max_spatial_level(),
                    other.min_spatial_level(),
                    other.max_spatial_index_cells(),
                    other.max_spatial_query_cells());
}

bool Indexer::Options::operator!=(const Options &other) const
{
    return !((*this) == other);
}

Indexer::Indexer(const Options &options)
{
    options_ = options;
    S2RegionTermIndexer::Options s2options;
    s2options.set_optimize_for_space(false); // Optimize for time.
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

vector<string> Indexer::terms(const TermStyle style, const S2Point &point)
{
    return generate_terms(style, point);
}

vector<string> Indexer::terms(const TermStyle style, const S2Region &region)
{
    return generate_terms(style, region);
}

vector<string> Indexer::terms(const TermStyle style, const Observation &observation)
{
    S2Polygon polygon;
    observation.as_polygon(polygon);
    return generate_terms(style, polygon);
}

vector<string> Indexer::terms(const TermStyle style, const Ephemeris &eph)
{
    return terms(style, eph, 0);
}

vector<string> Indexer::terms(const TermStyle style, const Ephemeris &eph, double padding)
{
    auto index = std::make_unique<MutableS2ShapeIndex>();
    for (auto &polygon : eph.as_polygons())
        index->Add(std::make_unique<S2Polygon::OwningShape>(std::move(polygon)));

    S2ShapeIndexBufferedRegion region(index.get(), S1ChordAngle::Degrees(padding / 60));
    return generate_terms(style, region);
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
    indexer_.mutable_options()->set_max_cells((style == index)
                                                  ? options_.max_spatial_index_cells()
                                                  : options_.max_spatial_query_cells());
    // spatial terms
    return (style == index)
               ? indexer_.GetIndexTerms(region, "")
               : indexer_.GetQueryTerms(region, "");
}