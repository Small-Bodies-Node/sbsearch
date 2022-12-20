#include "ephemeris.h"
#include "util.h"
#include "sbsearch_testing.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <stdexcept>
#include <vector>

#include <s2/s1angle.h>
#include <s2/s1chord_angle.h>
#include <s2/s2edge_distances.h>
#include <s2/s2latlng.h>
#include <s2/s2point.h>
#include <s2/s2polyline.h>
#include <s2/s2region_term_indexer.h>

using sbsearch::position_angle;
using std::cerr;
using std::cout;
using std::endl;
using std::make_unique;
using std::unique_ptr;
using std::vector;

namespace sbsearch
{
    Ephemeris::Ephemeris(const vector<S2Point> vertices, const vector<double> times)
    {
        assert(vertices.size() == times.size());

        for (auto time : times)
            times_.push_back(time);

        for (auto vertex : vertices)
            vertices_.push_back(vertex);

        num_vertices_ = vertices.size();
        num_segments_ = num_vertices_ - 1;
        isValid();
    }

    bool Ephemeris::isValid()
    {
        if (!testing::is_increasing(times_))
            throw std::runtime_error("Times must be monotonically increasing.");

        return true;
    }

    bool Ephemeris::is_equal(Ephemeris &other)
    {
        if (num_vertices() != other.num_vertices())
            return false;

        for (int i = 0; i < num_vertices(); i++)
        {
            if (vertex(i) != other.vertex(i))
                return false;

            if (time(i) != other.time(i))
                return false;
        }

        return true;
    }

    int Ephemeris::num_vertices()
    {
        return num_vertices_;
    }

    S2Point Ephemeris::vertex(const int k)
    {
        assert(k >= 0);
        assert(k < num_vertices());
        return vertices_[k];
    }

    vector<S2Point> Ephemeris::vertices()
    {
        vector<S2Point> v;
        for (int i = 0; i < num_vertices(); i++)
            v.push_back(vertices_[i]);
        return v;
    }

    double Ephemeris::time(const int k)
    {
        assert(k >= 0);
        assert(k < num_vertices());
        return times_[k];
    }

    vector<double> Ephemeris::times()
    {
        vector<double> t;
        for (int i = 0; i < num_vertices(); i++)
            t.push_back(time(i));
        return t;
    }

    int Ephemeris::num_segments()
    {
        return num_segments_;
    }

    Ephemeris Ephemeris::segment(const int k)
    {
        assert(k > -num_segments());
        assert(k < num_segments());

        const int offset = (k < 0) ? num_segments() - 1 : 0;
        vector<S2Point> v{vertices_[k + offset], vertices_[k + offset + 1]};
        vector<double> t{times_[k + offset], times_[k + offset + 1]};
        return Ephemeris(v, t);
    }

    vector<Ephemeris> Ephemeris::segments()
    {
        vector<Ephemeris> v;
        for (int i = 0; i < num_segments(); i++)
            v.push_back(segment(i));
        return v;
    }

    S2Polyline Ephemeris::as_polyline()
    {
        return S2Polyline(vertices());
    }

    S2Point Ephemeris::interpolate(const double mjd)
    {
        if ((mjd < *times_.begin()) | (mjd > *(times_.end() - 1)))
            throw std::runtime_error("Interpolation beyond ephemeris time range.");

        // find the nearest segment
        auto start = std::find_if(times_.begin(), times_.end(), [mjd](double t)
                                  { return (t > mjd); }) -
                     1;
        auto end = start + 1;

        // length of segment in time
        double dt = *end - *start;

        // interpolate to this fraction
        double frac = (mjd - *start) / dt;

        // this is the line we will interpolate
        S2Polyline segment(
            vector<S2Point>(vertices_.begin() + (start - times_.begin()),
                            vertices_.begin() + (end - times_.begin() + 1)));
        return segment.Interpolate(frac);
    }

    S2Point Ephemeris::extrapolate(const double distance, Ephemeris::Extrapolate direction)
    {
        S2Point a, b;
        if (direction == Ephemeris::Extrapolate::BACKWARDS)
        {
            a = vertices_[1];
            b = vertices_[0];
        }
        else
        {
            a = vertices_[num_vertices() - 2];
            b = vertices_[num_vertices() - 1];
        }
        double length = S1Angle(a, b).radians();
        S2Point extrapolated = S2::Interpolate(a, b, 1 + distance / length);
        return extrapolated.Normalize();
    }

    Ephemeris Ephemeris::subsample(const double mjd_start, const double mjd_stop)
    {
        vector<S2Point> vertices;
        vector<double> times;

        // start the line
        vertices.push_back(interpolate(mjd_start));
        times.push_back(mjd_start);

        // find any whole segments between start and end
        auto next = std::find_if(times_.begin(), times_.end(), [mjd_start](double t)
                                 { return (t > mjd_start); });
        auto last = std::find_if(times_.begin(), times_.end(), [mjd_stop](double t)
                                 { return (t > mjd_stop); }) -
                    1;

        if (*last > *next)
        {
            // there is at least one segment between start and end
            for (int i = (next - times_.begin()); i <= (last - times_.begin()); i++)
            {
                vertices.push_back(vertices_[i]);
                times.push_back(times_[i]);
            }
        }

        // end the line
        vertices.push_back(interpolate(mjd_stop));
        times.push_back(mjd_stop);

        return Ephemeris(vertices, times);
    }

    S2Polygon Ephemeris::pad(const vector<double> &para, const vector<double> &perp)
    {
        // para[1] through para[-2] are not used

        assert(para.size() == perp.size());
        assert(*std::max_element(para.begin(), para.end()) < 90 * DEG);
        assert(*std::max_element(perp.begin(), perp.end()) < 90 * DEG);

        vector<S2Point> spine;

        // extrapolate the vertices with para.
        spine.push_back(extrapolate(para[0], Extrapolate::BACKWARDS));
        for (auto vertex : vertices_)
            spine.push_back(vertex);
        spine.push_back(extrapolate(*(para.end() - 1), Extrapolate::FORWARDS));

        // extend perp to match spine, transform to angles
        vector<S1Angle> perp_angles;
        perp_angles.push_back(S1Angle::Radians(perp[0]));
        for (int i = 0; i < perp.size(); i++)
            perp_angles.push_back(S1Angle::Radians(perp[i]));
        perp_angles.push_back(S1Angle::Radians(*(perp.end() - 1)));

        // create pairs of parallel tracks with perp
        vector<S2Point> vertices;
        const int n = spine.size();
        for (int i = 0; i < n - 1; i++)
        {
            vertices.push_back(S2::GetPointToRight(spine[i], spine[i + 1], perp_angles[i]));
        }
        // calculate end using the reversed vector
        vertices.push_back(S2::GetPointToLeft(spine[n - 1], spine[n - 2], perp_angles[n - 1]));
        for (int i = n - 1; i > 0; i--)
        {
            vertices.push_back(S2::GetPointToRight(spine[i], spine[i - 1], perp_angles[i]));
        }
        vertices.push_back(S2::GetPointToLeft(spine[0], spine[1], perp_angles[n - 1]));

        S2Polygon polygon;
        makePolygon(vertices, polygon);
        return polygon;
    }

    S2Polygon Ephemeris::pad(const double para, const double perp)
    {
        vector<double> para_vector(num_vertices(), para);
        vector<double> perp_vector(num_vertices(), perp);
        return pad(para_vector, perp_vector);
    }

    // vector<string> Ephemeris::query_terms(S2RegionTermIndexer &indexer)
    // {
    //     return generate_terms(TermStyle::query, indexer);
    // }

    // vector<string> Ephemeris::query_terms(S2RegionTermIndexer &indexer, const vector<double> &para, vector<double> &perp)
    // {
    //     return generate_terms(TermStyle::query, indexer, para, perp);
    // }

    // vector<string> Ephemeris::query_terms(S2RegionTermIndexer &indexer, const double para, const double perp)
    // {
    //     vector<double> para_vector(num_vertices(), para);
    //     vector<double> perp_vector(num_vertices(), perp);
    //     return query_terms(indexer, para_vector, perp_vector);
    // }

    // vector<string> Ephemeris::index_terms(S2RegionTermIndexer &indexer)
    // {
    //     return generate_terms(TermStyle::index, indexer);
    // }

    // vector<string> Ephemeris::index_terms(S2RegionTermIndexer &indexer, const vector<double> &para, vector<double> &perp)
    // {
    //     return generate_terms(TermStyle::index, indexer, para, perp);
    // }

    // vector<string> Ephemeris::index_terms(S2RegionTermIndexer &indexer, const double para, const double perp)
    // {
    //     vector<double> para_vector(num_vertices(), para);
    //     vector<double> perp_vector(num_vertices(), perp);
    //     return index_terms(indexer, para_vector, perp_vector);
    // }

    // vector<string> Ephemeris::generate_terms(TermStyle style, S2RegionTermIndexer &indexer)
    // {
    //     vector<string> terms;
    //     for (auto eph : segments())
    //     {
    //         S2Polyline segment(eph.vertices());

    //         // Get query terms for the segment
    //         vector<string> spatial_terms;
    //         if (style == TermStyle::index)
    //             spatial_terms = indexer.GetIndexTerms(segment, "");
    //         else
    //             spatial_terms = indexer.GetQueryTerms(segment, "");

    //         // Get terms for the time
    //         vector<string> time_terms = mjd_to_time_terms(eph.time(0), eph.time(1));

    //         // Join query terms, each segment gets a time suffix
    //         for (auto time_term : time_terms)
    //             for (auto spatial_term : spatial_terms)
    //                 terms.push_back(spatial_term + "-" + time_term);
    //     }
    //     return terms;
    // }

    // vector<string> Ephemeris::generate_terms(TermStyle style, S2RegionTermIndexer &indexer, const vector<double> &para, vector<double> &perp)
    // {
    //     vector<string> terms;
    //     auto para_iterator = para.begin();
    //     auto perp_iterator = perp.begin();
    //     for (auto eph : segments())
    //     {
    //         unique_ptr<S2Polygon> polygon = eph.pad(vector<double>(para_iterator, para_iterator + 2), vector<double>(perp_iterator, perp_iterator + 2));

    //         // Get query terms for the segment
    //         vector<string> spatial_terms;
    //         if (style == TermStyle::index)
    //             spatial_terms = indexer.GetIndexTerms(*polygon, "");
    //         else
    //             spatial_terms = indexer.GetQueryTerms(*polygon, "");

    //         // Get terms for the time
    //         vector<string> time_terms = mjd_to_time_terms(eph.time(0), eph.time(1));

    //         // Join query terms, each segment gets a time suffix
    //         for (auto time_term : time_terms)
    //             for (auto spatial_term : spatial_terms)
    //                 terms.push_back(spatial_term + "-" + time_term);
    //     }
    //     return terms;
    // }
}