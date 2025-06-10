#include <optional>
#include <s2/s2cap.h>
#include <s2/s2debug.h>
#include <s2/s2polygon.h>
#include <s2/s2builderutil_s2polygon_layer.h>

#include "intersection.h"
#include "observation.h"

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

    bool intersects(const S2Polygon &polygon,
                    const S2Cap &cap,
                    const IntersectionType intersection_type)
    {
        bool result = false;
        switch (intersection_type)
        {
        case (ContainsPoint):
            result = polygon.Contains(cap.center());
            break;
        case (ContainsArea):
            result = (polygon.GetDistanceToBoundary(cap.center()) > cap.radius().ToAngle() & polygon.Contains(cap.center()));
            break;
        case (IntersectsArea):
            result = polygon.GetDistance(cap.center()) < cap.radius().ToAngle();
            break;
        case (ContainedByArea):
            // only testing loop[0]; sbsearch does not use multiple loops
            const S2Loop *loop = polygon.loop(0);

            // check that each vertex is contained; immediately end loop if any
            // vertex is not
            result = true;
            for (int i = 0; i < loop->num_vertices(); i++)
            {
                result = cap.InteriorContains(loop->vertex(i));
                if (!result)
                    break;
            }
            break;
        }
        return result;
    }

    bool intersects(const S2Polygon &polygon1,
                    const S2Polygon &polygon2,
                    const IntersectionType intersection_type)
    {
        bool result = false;

        switch (intersection_type)
        {
        case (ContainsCenter):
            result = polygon1.Contains(polygon2.GetCentroid().Normalize());
            break;
        case (ContainsArea):
            result = polygon1.Contains(polygon2);
            break;
        case (IntersectsArea):
            result = polygon1.Intersects(polygon2);
            break;
        case (ContainedByArea):
            result = polygon2.Contains(polygon1);
            break;
        }
        return result;
    }

    bool intersects(const Observation &observation,
                    const optional<double> mjd_start,
                    const optional<double> mjd_stop)
    {
        if (mjd_start.has_value() & ((observation.mjd_stop() < mjd_start)))
            return false;

        if (mjd_stop.has_value() & ((observation.mjd_start() > mjd_stop)))
            return false;

        return true;
    }

    bool contains(const Observation &observation,
                  const S2Point &point,
                  const optional<double> mjd_start,
                  const optional<double> mjd_stop)
    {
        if (!intersects(observation, mjd_start, mjd_stop))
            return false;

        S2Polygon fov;
        fov.set_s2debug_override(S2Debug::DISABLE);
        observation.as_polygon(fov);
        return fov.Contains(point);
    }

    bool intersects(const Observation &observation,
                    const S2Cap &cap,
                    const IntersectionType intersection_type,
                    const optional<double> mjd_start,
                    const optional<double> mjd_stop)
    {
        if (!intersects(observation, mjd_start, mjd_stop))
            return false;

        S2Polygon fov;
        fov.set_s2debug_override(S2Debug::DISABLE);
        observation.as_polygon(fov);
        return intersects(fov, cap, intersection_type);
    }

    bool intersects(const Observation &observation,
                    const S2Polygon &area,
                    const IntersectionType intersection_type,
                    const optional<double> mjd_start,
                    const optional<double> mjd_stop)
    {
        if (!intersects(observation, mjd_start, mjd_stop))
            return false;

        S2Polygon fov;
        fov.set_s2debug_override(S2Debug::DISABLE);
        observation.as_polygon(fov);
        return intersects(fov, area, intersection_type);
    }

    bool intersects(const Observation &observation,
                    const Ephemeris &ephemeris,
                    const optional<double> padding,
                    const optional<double> mjd_start,
                    const optional<double> mjd_stop)
    {
        if (!intersects(observation, mjd_start, mjd_stop))
            return false;

        S2Polygon fov;
        fov.set_s2debug_override(S2Debug::DISABLE);
        observation.as_polygon(fov);

        if (ephemeris.options().use_uncertainty)
        {
            auto polygons = ephemeris.as_polygons(padding.value_or(0));
            for (auto const &polygon : polygons)
            {
                if (intersects(fov, *polygon, IntersectsArea))
                    return true;
            }
        }
        else
        {
            auto line = ephemeris.subsample(
                                     std::max(mjd_start.value_or(0), observation.mjd_start()),
                                     std::min(mjd_stop.value_or(100000), observation.mjd_stop()))
                            .as_polyline();
            return fov.Intersects(line);
        }
        return false;
    }
}
