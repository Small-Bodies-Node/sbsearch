# Licensed with the 3-clause BSD license.  See LICENSE for details.
"""Interior testing for spherical geometry."""
from libc.math cimport sin, cos, acos, fabs, M_PI
import numpy as np


def spherical_distance(double ra1, double dec1,
                       double ra2, double dec2):
    return acos(sin(dec1) * sin(dec2) + cos(dec1) * cos(dec2) * cos(fabs(ra1 - ra2)))


def spherical_triangle_area(double a, double b, double c):
    """Area of spherical triangle."""
    cdef double ca, cb, cc, sa, sb, sc, A, B, C
    ca = cos(a)
    cb = cos(b)
    cc = cos(c)
    sa = sin(a)
    sb = sin(b)
    sc = sin(c)
    A = acos(max(-1, min(1, (ca - cb * cc) / (sb * sc))))
    B = acos(max(-1, min(1, (cb - cc * ca) / (sc * sa))))
    C = acos(max(-1, min(1, (cc - ca * cb) / (sa * sb))))

    return A + B + C - M_PI


def interior_test(point, corners, rtol=1e-5):
    """Test if point is interior to corners assuming spherical geometry.


    Parameters
    ----------
    point : RADec
        Point to test.

    corners : RADec
        Points describing a spherical rectangle.

    rtol : float, optional
        Relative tolerance for area comparison.


    Returns
    -------
    interior :
        ``True`` if the point falls inside the rectangle.

    """

    cdef double d[4], segments[5], area_r, area_p
    cdef int ii, i, j, k

    for ii in range(3):
        d[ii] = spherical_distance(corners.ra[0], corners.dec[0],
                                   corners.ra[ii + 1], corners.dec[ii + 1])

    # 0, k and i, j are opposite corners
    i, j, k = np.argsort((d[0], d[1], d[2])) + 1

    segments[0] = spherical_distance(corners.ra[0], corners.dec[0],
                                     corners.ra[i], corners.dec[i])
    segments[1] = spherical_distance(corners.ra[i], corners.dec[i],
                                     corners.ra[k], corners.dec[k])
    segments[2] = spherical_distance(corners.ra[k], corners.dec[k],
                                     corners.ra[j], corners.dec[j])
    segments[3] = spherical_distance(corners.ra[j], corners.dec[j],
                                     corners.ra[0], corners.dec[0])
    segments[4] = spherical_distance(corners.ra[j], corners.dec[j],
                                     corners.ra[i], corners.dec[i])

    # measure area of the rectangle
    area_r = (spherical_triangle_area(segments[0], segments[1], segments[4])
              + spherical_triangle_area(segments[2], segments[3], segments[4]))

    # measure area of all triangles made by point and rectangle segments
    for ii in range(4):
        d[ii] = spherical_distance(point.ra, point.dec,
                                   corners.ra[ii], corners.dec[ii])
    area_p = (spherical_triangle_area(segments[0], d[0], d[i])
              + spherical_triangle_area(segments[1], d[i], d[k])
              + spherical_triangle_area(segments[2], d[k], d[j])
              + spherical_triangle_area(segments[3], d[j], d[0]))

    return np.isclose(area_r, area_p, rtol=rtol)
