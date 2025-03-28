#ifndef INTERSECTION_H_
#define INTERSECTION_H_

#include <s2/s2cap.h>
#include <s2/s2polygon.h>

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
    bool intersects(const S2Polygon &polygon, const S2Cap &area, const IntersectionType intersection_type);

    // Test for intersection between two polygons.
    bool intersects(const S2Polygon &polygon, const S2Polygon &area, const IntersectionType intersection_type);
}

#endif // INTERSECTION_H_
