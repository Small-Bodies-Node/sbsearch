#ifndef INTERSECTION_H_
#define INTERSECTION_H_

#include <optional>
#include <s2/s2cap.h>
#include <s2/s2polygon.h>
#include "observation.h"
#include "ephemeris.h"

namespace sbsearch
{
    // Intersection types
    enum IntersectionType
    {
        ContainsPoint,
        ContainsArea,
        IntersectsArea,
        ContainedByArea,
        ContainsCenter = ContainsPoint
    };

    std::istream &operator>>(std::istream &in, IntersectionType &intersection_type);

    // Test for intersection between a polygon and a spherical cap.
    bool intersects(const S2Polygon &polygon,
                    const S2Cap &cap,
                    const IntersectionType intersection_type);

    // Test for intersection between two polygons.
    bool intersects(const S2Polygon &polygon1,
                    const S2Polygon &polygon2,
                    const IntersectionType intersection_type);

    // Test if observation intersects given time range
    bool intersects(const Observation &observation,
                    const optional<double> mjd_start,
                    const optional<double> mjd_stop);

    // Test if observation contains point, with optional time limits.
    bool contains(const Observation &observation,
                  const S2Point &point,
                  const optional<double> mjd_start = std::nullopt,
                  const optional<double> mjd_stop = std::nullopt);

    // Test for intersection between observation and spherical cap, with optional time limits.
    bool intersects(const Observation &observation,
                    const S2Cap &cap,
                    const IntersectionType intersection_type,
                    const optional<double> mjd_start = std::nullopt,
                    const optional<double> mjd_stop = std::nullopt);

    // Test for intersection between observation and polygon, with optional time limits.
    bool intersects(const Observation &observation,
                    const S2Polygon &area,
                    const IntersectionType intersection_type,
                    const optional<double> mjd_start = std::nullopt,
                    const optional<double> mjd_stop = std::nullopt);

    // Test for intersection between observation and ephemeris, with optional
    // padding and time limits.
    bool intersects(const Observation &observation,
                    const Ephemeris &ephemeris,
                    const optional<double> padding = std::nullopt,
                    const optional<double> mjd_start = std::nullopt,
                    const optional<double> mjd_stop = std::nullopt);
}

#endif // INTERSECTION_H_
