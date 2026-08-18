/** S2Polygon generation functions for sbsearch objects.
 *
 * These functions are generally designed for speed.  Creating/destroying many
 * S2Polygon objects seems to be slow.
 *
 */
#ifndef SBS_POLYGONS_H_
#define SBS_POLYGONS_H_

#include <vector>
#include <s2/s2polygon.h>

#include "ephemeris/ephemeris.h"
#include "observation.h"

using sbsearch::ephemeris::Ephemeris;
using std::vector;

namespace sbsearch
{
    // Make a single simple polygon (no edge crossings, etc.).  The orientation
    // of the vertices will forced to make a polygon of area smaller than 2 pi.
    void make_polygon_simple(const vector<S2Point> &vertices, S2Polygon &polygon);

    // Make a polygon from a set of vertices will all edge and error checks.
    void make_polygon(const vector<S2Point> &vertices, S2Polygon &polygon);

    // Pad around the polygon by padding in arcmin.  Padding must be
    // > 0 or else the original polygon will be returned unmodified.
    void padded_polygon(const S2Polygon &polygon,
                        const double padding,
                        S2Polygon &result);

    /**
    Convert the ephemeris into a vector of polygons, with optional
    padding in arcmin.

    The area will depend on the `use_uncertainty` option, but the padding
    around the ephemeris will be at least 2", based on results in the testing suite.

    The ephemeris is described by connecting parallelograms that
    circumscribe ellipses with semi-major axes `a`, semi-minor axes `b`,
    with `a` aligned along angle `theta` (E of N).
    */
    vector<std::unique_ptr<S2Polygon>> make_polygons(const Ephemeris &eph,
                                                     const bool use_uncertainty,
                                                     double padding = 0);

    // Create a polygon from this observation's field of view, with optional
    // validation checks.
    void make_polygon(const Observation &observation,
                      S2Polygon &polygon,
                      const bool verify = false);
}

#endif
