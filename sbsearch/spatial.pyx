# Licensed with the 3-clause BSD license.  See LICENSE for details.
# cython: language_level=3
# distutils: language = c++

import cython
from libcpp.memory cimport unique_ptr, make_unique
from libcpp.utility cimport move
from libc.math cimport atan2, asin, sin, cos, tan, M_PI, M_PI_2
from libcpp.vector cimport vector

import numpy as np
import pywraps2 as s2
cimport numpy as np
from .spatial cimport (
    kAvgEdge,
    _build_polygon,
    _verify_build_polygon,
    _polygon_contains_point,
    _polygon_intersects_line,
    _polygon_intersects_about_line,
    _polygon_intersects_polygon,
    _polygon_intersects_cap,
    _position_angle,
    _offset_by,
)

cdef S2Point NORTH_CELESTIAL_POLE = S2Point(0, 0, 1)
cdef S2Point VERNAL_EQUINOX = S2Point(1, 0, 0)


class PolygonBuildError(Exception):
    pass


def closest_level(edge_length):
    """edge_length in radias"""
    return kAvgEdge.GetClosestLevel(edge_length)


def position_angle(double ra1, double dec1, double ra2, double dec2):
    """ra, dec in radians"""
    return _position_angle(ra1, dec1, ra2, dec2)


def offset_by(double ra, double dec, double pa, double rho):
    """ra, dec, pa, rho in radians"""
    cdef double new_ra = 0, new_dec = 0
    _offset_by(ra, dec, pa, rho, new_ra, new_dec)
    return new_ra, new_dec


def verify_build_polygon(ra, dec):
    """ra dec in radians, polygon is not closed"""
    cdef double[::1] ra_memview = ra.copy(order="C")
    cdef double[::1] dec_memview = dec.copy(order="C")
    try:
        return _verify_build_polygon(&ra_memview[0], &dec_memview[0], ra_memview.shape[0])
    except ValueError as exc:
        raise PolygonBuildError(str(exc))


def polygon_intersects_line(poly_ra, poly_dec, line_ra, line_dec,
                            double line_start=0, double line_stop=1):
    """Test if the polygon intersects line.
    

    Parameters
    ----------
    poly_ra, poly_dec : ndarray
        Polygon RA and Dec, radians.

    line_ra, line_dec : ndarray
        Line RA and Dec, radians.

    line_start, line_stop : float
        Test a sub-portion of the line, indicated by the given line fractions
        (linear interpolation).

    """

    cdef double[::1] poly_ra_memview = poly_ra.copy(order="C")
    cdef double[::1] poly_dec_memview = poly_dec.copy(order="C")
    cdef double[::1] line_ra_memview = line_ra.copy(order="C")
    cdef double[::1] line_dec_memview = line_dec.copy(order="C")

    return _polygon_intersects_line(&poly_ra_memview[0], &poly_dec_memview[0], poly_ra_memview.shape[0],
                                    &line_ra_memview[0], &line_dec_memview[0], line_ra_memview.shape[0],
                                    line_start, line_stop)


def polygon_intersects_about_line(poly_ra, poly_dec, line_ra, line_dec, a, b,                                  
                                  double line_start=0, double line_stop=1):
    """Test for intersection between polygon and region about a line.


    Parameters
    ----------
    poly_ra, poly_dec : ndarray
        Polygon RA and Dec, radians.

    line_ra, line_dec : ndarray
        Line RA and Dec, radians.

    a : ndarray
        Padding parallel to the line segements (only the first and last values are used), radians.

    b : ndarray
        Padding perpendicular to the line segements, radians.

    line_start, line_stop : float
        Test a sub-portion of the line, indicated by the given line fractions
        (linear interpolation).




    Returns
    -------
    intersects : bool
  
    """

    cdef double[::1] poly_ra_memview = poly_ra.copy(order="C")
    cdef double[::1] poly_dec_memview = poly_dec.copy(order="C")
    cdef double[::1] line_ra_memview = line_ra.copy(order="C")
    cdef double[::1] line_dec_memview = line_dec.copy(order="C")
    cdef double[::1] a_memview = a.copy(order="C")
    cdef double[::1] b_memview = b.copy(order="C")

    return _polygon_intersects_about_line(
        &poly_ra_memview[0], &poly_dec_memview[0], poly_ra_memview.shape[0],
        &line_ra_memview[0], &line_dec_memview[0], line_ra_memview.shape[0],
        &a_memview[0], &b_memview[0], line_start, line_stop)


def polygon_intersects_polygon(ra1, dec1, ra2, dec2):
    """Test intersection between two polygons.
    
    
    Parameters
    ----------
    ra1, dec1, ra2, dec2 : array
        Vertices of the polygons.  The polygons will be closed.  Units of radians.

    """

    cdef double[::1] ra1_memview = ra1.copy(order="C")
    cdef double[::1] dec1_memview = dec1.copy(order="C")
    cdef double[::1] ra2_memview = ra2.copy(order="C")
    cdef double[::1] dec2_memview = dec2.copy(order="C")

    return _polygon_intersects_polygon(&ra1_memview[0], &dec1_memview[0], ra1_memview.shape[0],
                                       &ra2_memview[0], &dec2_memview[0], ra2_memview.shape[0])


def polygon_contains_point(poly_ra, poly_dec, double point_ra, double point_dec):
    """Test if polygon contains a point.
    
    
    Parameters
    ----------
    poly_ra, poly_dec : ndarray
        Vertices of the polygon.  The polygon will be closed.  Units of radians.

    point_ra, point_dec : float
        The point to test, radians.

    """

    cdef double[::1] poly_ra_memview = poly_ra.copy(order="C")
    cdef double[::1] poly_dec_memview = poly_dec.copy(order="C")

    return _polygon_contains_point(&poly_ra_memview[0], &poly_dec_memview[0], poly_ra_memview.shape[0],
                                   point_ra, point_dec)


def polygon_intersects_cap(poly_ra, poly_dec, double point_ra, double point_dec,
                           double radius, int intersection_type):
    """Test if polygon intersects a cap.
    
    
    Parameters
    ----------
    poly_ra, poly_dec : ndarray
        Vertices of the polygon.  The polygon will be closed.  Units of radians.

    point_ra, point_dec : float
        The cap center, radians.

    radius : float
        The cap radius, radians.

    intersection_type : int
        Intersection type, see libspatial.cpp.

    """

    cdef double[::1] poly_ra_memview = poly_ra.copy(order="C")
    cdef double[::1] poly_dec_memview = poly_dec.copy(order="C")

    return _polygon_intersects_cap(&poly_ra_memview[0], &poly_dec_memview[0], poly_ra_memview.shape[0],
                                   point_ra, point_dec, radius, intersection_type)


def term_to_cell_vertices(term):
    """Convert spatial index term to S2 cell vertices.
    

    Parameters
    ----------
    term : string


    Returns
    -------
    ra, dec : ndarray
        Radians.
        
    """

    cell = s2.s2cell(s2.S2CellId.FromToken(term))

    ra, dec = np.empty((2, 4))
    for i in range(4):
        latlng = s2.S2LatLng(cell.GetVertex(i)).Normalized()
        ra[i] = latlng.lng().radians()
        dec[i] = latlng.lat().radians()

    return ra, dec


class SpatialIndexer:
    """Spatial indexer.

    Parameters
    ----------
    min_edge_length : double
    max_edge_length : double
        Minimum and maximum edge length to index, radians.

    max_cells : int, optional
        Maximum number of cells generated when approximating each region
        (more cells may be generated depending on min/max level).

    """

    def __init__(self, min_edge_length, max_edge_length, max_cells: int = 8):
        self._indexer = s2.S2RegionTermIndexer()
        self._indexer.set_max_cells(max_cells)
        self._indexer.set_max_level(kAvgEdge.GetClosestLevel(min_edge_length))
        self._indexer.set_min_level(kAvgEdge.GetClosestLevel(max_edge_length))

    @staticmethod
    def vertices_to_loop(ra, dec):
        loop = s2.S2Loop()
        loop.Init([s2.S2LatLng.FromRadians(*v).Normalized().ToPoint() for v in zip(dec, ra)])
        loop.Normalize()
        return loop

    def index_polygon(self, ra, dec):
        """Index the given polygon.


        Parameters
        ----------
        ra, dec : ndarray
            Polygon vertices, radians.  Order  (CW vs CCW) does not matter.  It is assumed the
            polygon should be smaller than 1/2 the sphere.


        Returns
        -------
        terms : list of strings


        """

        n = ra.shape[0]
        if dec.shape[0] != n:
            raise ValueError('ra and dec have different lengths')

        polygon = s2.S2Polygon(self.vertices_to_loop(ra, dec))
 
        terms = self._indexer.GetIndexTerms(polygon, "")
        return terms

    def query_point(self, double ra, double dec):
        """Query terms for a point.
        
        
        Parameters
        ----------
        ra, dec : ndarray
            The point to query, radians.


        Returns
        -------
        terms : list of strings
        
        """
        
        p = s2.S2LatLng.FromRadians(dec, ra).Normalized().ToPoint()
        terms = self._indexer.GetQueryTerms(p, "")
        return terms

    def query_cap(self, double ra, double dec, double radius):
        """Query terms for a spherical cap.
        

        Parameters
        ----------
        ra, dec : float
            The point to query, radians.

        radius : float
            The radius of the cap, radians.


        Returns
        -------
        terms : list of strings

        """

        center = s2.S2LatLng.FromRadians(dec, ra).Normalized().ToPoint()
        cap = s2.S2Cap(center, s2.S1Angle.Radians(radius))
        terms = self._indexer.GetQueryTerms(cap, "")
        return terms

    @staticmethod
    def vertices_to_polyline(ra, dec):
        vertices = [s2.S2LatLng.FromRadians(*v).Normalized() for v in zip(dec, ra)]
        line = s2.S2Polyline()
        line.InitFromS2LatLngs(vertices)
        return line

    def query_line(self, ra, dec):
        """Query terms for a region defined by a multi-segmented line.


        Parameters
        ----------
        ra, dec : ndarray
            Coordinates defining line vertices, radians.


        Returns
        -------
        terms : list of strings

        """

        n = ra.shape[0]
        if dec.shape[0] != n:
            raise ValueError('ra and dec have different lengths')

        line = self.vertices_to_polyline(ra, dec)

        terms = self._indexer.GetQueryTerms(line, "")
        return terms

    def query_polygon(self, ra, dec):
        """Query terms for a polygon.


        Parameters
        ----------
        ra, dec : ndarray
            Polygon vertices, radians.  Order does not matter.  It is assumed
            the polygon should be smaller than 1/2 the sphere.


        Returns
        -------
        terms : list of strings

        """

        n = ra.shape[0]
        if dec.shape[0] != n:
            raise ValueError('ra and dec have different lengths')

        polygon = s2.S2Polygon(self.vertices_to_loop(ra, dec))

        terms = self._indexer.GetQueryTerms(polygon, "")
        return terms

    def query_about_line(self, ra, dec, a, b):
        """Query terms for a region surrounding a multi-segmented line.


        Parameters
        ----------
        ra, dec : ndarray
            Coordinates defining line vertices, radians.

        a, b : ndarray
            Angular size of the region: ``a`` is along the line, ``b`` 
            is perpendicular to it, radians.  Only the first and last
            ``a`` are considered.


        Returns
        -------
        terms : list of byte strings
            Query terms for the polygon.

        poly_ra, poly_dec : np.ndarray
            The generated polygon vertices.

        """

        n = ra.shape[0]
        n_dec = dec.shape[0]
        n_a = a.shape[0]
        n_b = b.shape[0]
        if any([n != n_dec, n_dec != n_a, n_a != n_b]):
            raise ValueError('All arrays must have equal length')

        # transform the line into a polygon
        # the spine is the input line extended by + / -a
        # the polygon is the region betweeh spine + b and spine - b
        pa = []
        for i in range(n - 1):
            pa.append(position_angle(ra[i], dec[i], ra[i + 1], dec[i + 1]))
        pa.insert(0, pa[0])
        pa.append(pa[-1])
        pa.append(pa[-1])

        spine = []
        for i in range(n):
            spine.append((ra[i], dec[i]))
        spine.insert(0, offset_by(ra[0], dec[0], pa[0], -a[0]))
        spine.append(offset_by(ra[-1], dec[-1], pa[-1], a[-1]))

        # Construct query region polygon. In S2, the interior of a line is on
        # the left.  This is CCW from each edge for small loops.  S2 is designed
        # for the surface of the Earth, but equatorial coordinates on the sky
        # are a mirror of the Earth's coordinates.  Therefore, the interior is
        # on the right (the CW edge).  Thus, we first generate the edges along
        # PA + pi / 2, then return to the first vertex along PA - pi / 2.

        vertices = []
        for i in range(n + 2):
            j = min(max(0, i - 1), n - 1)
            vertices.append(offset_by(spine[i][0], spine[i][1], pa[i] + M_PI_2, b[j]))
            vertices.insert(0, offset_by(spine[i][0], spine[i][1], pa[i] - M_PI_2, b[j]))

        poly_ra, poly_dec = np.array(vertices).T
        return self.query_polygon(poly_ra, poly_dec), poly_ra, poly_dec