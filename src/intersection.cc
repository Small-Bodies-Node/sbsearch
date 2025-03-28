#include <s2/s2cap.h>
#include <s2/s2polygon.h>

#include "intersection.h"

namespace sbsearch
{
    std::istream &operator>>(std::istream &in, sbsearch::IntersectionType &intersection_type)
    {
        std::string token;
        in >> token;
        std::transform(token.begin(), token.end(), token.begin(),
                       [](unsigned char c)
                       { return std::tolower(c); });

        if ((token == "containspoint") || (token == "containscenter"))
            intersection_type = IntersectionType::ContainsPoint;
        else if (token == "containsarea")
            intersection_type = IntersectionType::ContainsArea;
        else if (token == "intersectsarea")
            intersection_type = IntersectionType::IntersectsArea;
        else if (token == "containedbyarea")
            intersection_type = IntersectionType::ContainedByArea;
        else
            in.setstate(std::ios_base::failbit);
        return in;
    }

    bool intersects(const S2Polygon &polygon, const S2Cap &area, const IntersectionType intersection_type)
    {
        bool result = false;
        switch (intersection_type)
        {
        case (ContainsPoint):
            result = polygon.Contains(area.center());
            break;
        case (ContainsArea):
            result = (polygon.GetDistanceToBoundary(area.center()) > area.radius().ToAngle() & polygon.Contains(area.center()));
            break;
        case (IntersectsArea):
            result = polygon.GetDistance(area.center()) < area.radius().ToAngle();
            break;
        case (ContainedByArea):
            // only testing loop[0]; sbsearch does not use multiple loops
            const S2Loop *loop = polygon.loop(0);

            // check that each vertex is contained; immediately end loop if any
            // vertex is not
            result = true;
            for (int i = 0; i < loop->num_vertices(); i++)
            {
                result = area.InteriorContains(loop->vertex(i));
                if (!result)
                    break;
            }
            break;
        }
        return result;
    }

    bool intersects(const S2Polygon &polygon, const S2Polygon &area, const IntersectionType intersection_type)
    {
        bool result = false;

        switch (intersection_type)
        {
        case (ContainsCenter):
            result = polygon.Contains(area.GetCentroid().Normalize());
            break;
        case (ContainsArea):
            result = polygon.Contains(area);
            break;
        case (IntersectsArea):
            result = polygon.Intersects(area);
            break;
        case (ContainedByArea):
            result = area.Contains(polygon);
            break;
        }
        return result;
    }
}
