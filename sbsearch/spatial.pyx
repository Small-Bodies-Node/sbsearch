# Licensed with the 3-clause BSD license.  See LICENSE for details.
# cython: language_level=3

import cython
from libcpp.memory cimport unique_ptr, make_unique
from libcpp.utility cimport move
from libc.math cimport atan2, asin, sin, cos, tan, M_PI, M_PI_2
from libcpp.vector cimport vector
import numpy as np
cimport numpy as np
from .spatial cimport (S2RegionTermIndexer, kAvgEdge, S2LatLng,
    S2LatLngRect, S2Polyline, S2Point, S2Loop, S2Polygon)

cdef S2Point NORTH_CELESTIAL_POLE = S2Point(0, 0, 1)
cdef S2Point VERNAL_EQUINOX = S2Point(1, 0, 0)


def closest_level(min_edge_length):
    return kAvgEdge.GetClosestLevel(min_edge_length)

def position_angle(ra1, dec1, ra2, dec2):
    return _position_angle(ra1, dec1, ra2, dec2)


cdef _position_angle(double ra1, double dec1, double ra2, double dec2):
    """Meeus 1998, Astronomical Algorithms, p116"""
    dra = ra2 - ra1
    return atan2(sin(dra), cos(dec1) * tan(dec2) - sin(dec1) * cos(dra))


def offset_by(ra, dec, pa, rho):
    return _offset_by(ra, dec, pa, rho)


cdef _offset_by(double ra, double dec, double pa, double rho):
    """Algorithm based on astropy.coordinates.angle_utilities.offset_by."""

    cdef cos_a = np.cos(rho)
    cdef sin_a = np.sin(rho)
    # note: c is measured from the pole, dec from equator
    cdef cos_c = np.sin(dec)
    cdef sin_c = np.cos(dec)
    cdef cos_B = np.cos(pa)
    cdef sin_B = np.sin(pa)

    cdef cos_b = cos_c * cos_a + sin_c * sin_a * cos_B
    cdef xsin_A = sin_a * sin_B * sin_c
    cdef xcos_A = cos_a - cos_b * cos_c

    cdef A = np.arctan2(xsin_A, xcos_A)

    if sin_c < 1e-12:
        # handle the pole as a very small angle
        A = M_PI_2 + cos_c * (M_PI_2 - pa)

    return ra + A, asin(cos_b)

def polygon_intersects_line(double[:] poly_ra, double[:] poly_dec,
                            double[:] line_ra, double[:] line_dec,
                            double line_start=0, double line_stop=1):
    """Test if the polygon intersects the line.
    

    Parameters
    ----------
    poly_ra, poly_dec : ndarray
        Polygon RA and Dec, radians.

    line_ra, line_dec : ndarray
        Line RA and Dec, radians.

    line_start, line_stop : float
        Test a sub-portion of the line, indicated by the given line fractions
        (linear interpolation).


    Returns
    -------
    intersects : bool
    
    """

    cdef int n = poly_ra.shape[0]
    if poly_dec.shape[0] != n:
        raise ValueError('ra and dec have different lengths')

    cdef vector[S2Point] poly_vertices
    cdef int i
    for i in range(n):
        poly_vertices.push_back(
            S2LatLng.FromRadians(poly_dec[i], poly_ra[i])
            .Normalized().ToPoint()
        )

    cdef unique_ptr[S2Loop] loop
    loop = make_unique[S2Loop](poly_vertices)
    loop.get().Normalize()
    cdef S2Polygon poly
    poly.Init(move(loop))

    n = line_ra.shape[0]
    if line_dec.shape[0] != n:
        raise ValueError('ra and dec have different lengths')

    cdef vector[S2Point] line_vertices
    for i in range(n):
        line_vertices.push_back(
            S2LatLng.FromRadians(line_dec[i], line_ra[i])
            .Normalized().ToPoint()
        )
    cdef S2Polyline line = S2Polyline(line_vertices)

    # interpolate to a sub-set?
    cdef S2Point start = line.Interpolate(line_start)
    cdef S2Point stop = line.Interpolate(line_stop)
    cdef int next_vertex, last_vertex
    cdef vector[S2Point] new_line_vertices
    if line_start != 0 or line_stop != 1:
        if line_stop <= line_start:
            raise ValueError('line_stop <= line_start')

        line.GetSuffix(line_start, &next_vertex)
        line.GetSuffix(line_stop, &last_vertex)
        # build the new line
        new_line_vertices.push_back(start)
        if next_vertex != last_vertex:
            for i in range(last_vertex, next_vertex):
                new_line_vertices.push_back(line_vertices[i])
        new_line_vertices.push_back(stop)

        line = S2Polyline(new_line_vertices)
   
    return poly.Intersects(line)

def polygon_string_intersects_line(s, *args, **kwargs):
    """Test if the polygon intersects line.
    

    Parameters
    ----------
    s : str
        Comma-separated RA:Dec pairs in degrees, e.g., "1:1, 1:2, 2:1".

    *args, **kwargs
        See polygon_intersects_line
    
    """

    coords = np.radians(np.array([c.split(':') for c in s.split(',')], float))
    return polygon_intersects_line(coords[:, 0], coords[:, 1], *args, **kwargs)


def polygon_intersects_about_line(double[:] poly_ra, double[:] poly_dec,
                                  double[:] line_ra, double[:] line_dec,
                                  double[:] a, double[:] b,
                                  double line_start=0, double line_stop=1):
    """Test for intersection between polygon and region about a line.
    

    Parameters
    ----------
    poly_ra, poly_dec : ndarray
        Polygon RA and Dec, radians.

    line_ra, line_dec : ndarray
        Line RA and Dec, radians.

    line_start, line_stop : float
        Test a sub-portion of the line, indicated by the given line fractions.

    a, b : ndarray
        Angular size of the region: ``a`` is along the line, ``b`` 
        is perpendicular to it, radians.  Only the first and last
        ``a`` are considered.

    line_start, line_stop : float
        Test a sub-portion of the line, indicated by the given line fractions
        (linear interpolation).


    Returns
    -------
    intersects : bool

    """

    cdef int i, j

    # validate line inputs; polygon is validated in polygon_intersects_polygon
    cdef int n = line_ra.shape[0]
    if any([n != line_dec.shape[0], n != a.shape[0], n != b.shape[0]]):
        raise ValueError('Line arrays must all have equal length')

    # define the line or sub-segment thereof
    cdef vector[S2Point] line_vertices
    for i in range(n):
        line_vertices.push_back(
            S2LatLng.FromRadians(line_dec[i], line_ra[i])
            .Normalized().ToPoint()
        )
    cdef S2Polyline line = S2Polyline(line_vertices)

    # interpolate to a sub-segment?
    cdef S2Point start = line.Interpolate(line_start)
    cdef S2Point stop = line.Interpolate(line_stop)
    cdef int next_vertex, last_vertex
    cdef vector[S2Point] new_line_vertices
    if line_start != 0 or line_stop != 1:
        if line_stop <= line_start:
            raise ValueError('line_stop <= line_start')

        line.GetSuffix(line_start, &next_vertex)
        line.GetSuffix(line_stop, &last_vertex)
        # build the new line
        new_line_vertices.push_back(start)
        if next_vertex != last_vertex:
            for i in range(last_vertex, next_vertex):
                new_line_vertices.push_back(line_vertices[i])
        new_line_vertices.push_back(stop)

        line = S2Polyline(new_line_vertices)

    # transform the line into a polygon

    # the spine is the input line extended by +/-a
    # the polygon is the region betweeh spine+b and spine-b
    cdef double[:,:] spine = np.empty((n + 2, 2))
    cdef double[:] _ra = np.empty(2 * (n + 2))
    cdef double[:] _dec = np.empty(2 * (n + 2))

    # direction of each vector as position angle, including extension
    cdef double[:] pa = np.empty(n + 2)
    for i in range(n - 1):
        pa[i + 1] = _position_angle(line_ra[i], line_dec[i], line_ra[i + 1], line_dec[i + 1])
    pa[0] = pa[1]
    pa[n] = pa[n - 1]
    pa[n + 1] = pa[n]

    spine[0, 0], spine[0, 1] = _offset_by(line_ra[0], line_dec[0], pa[0], -a[0])
    for i in range(n):
        spine[i + 1, 0] = line_ra[i]
        spine[i + 1, 1] = line_dec[i]
    spine[n + 1, 0], spine[n + 1, 1] = _offset_by(
        line_ra[n - 1], line_dec[n - 1], pa[n - 1], a[n - 1])

    # pad out b to match the number of spine vertices
    cdef double[:] _b = np.r_[b[0], b, b[n - 1]]

    # construct polygon
    # first half of the region 
    for i in range(n + 2):
        _ra[i], _dec[i] = _offset_by(spine[i, 0], spine[i, 1], pa[i] + M_PI_2, _b[i])

    # second half
    for i in range(n + 2):
        j = n + 1 - i
        _ra[n + 2 + i], _dec[n + 2 + i] = _offset_by(spine[j, 0], spine[j, 1], pa[j] - M_PI_2, _b[j])

    return polygon_intersects_polygon(poly_ra, poly_dec, _ra, _dec)


def polygon_string_intersects_about_line(s, *args, **kwargs):
    """Test for intersection between polygon and region about a line.


    Parameters
    ----------
    s : str
        Comma-separated RA:Dec pairs in degrees, e.g., "1:1, 1:2, 2:1".

    *args, **kwargs
        Line arguments passed to ``polygon_intersects_about_line``.


    Returns
    -------
    intersects : bool
  
    """
    coords = np.radians(np.array([c.split(':') for c in s.split(',')], float))
    return polygon_intersects_about_line(coords[:, 0], coords[:, 1], *args, **kwargs)


def polygon_intersects_polygon(double[:] ra1, double[:] dec1,
                               double[:] ra2, double[:] dec2):
    """Test for polygon intersection."""
    
    cdef int n = ra1.shape[0]
    if dec1.shape[0] != n:
        raise ValueError('ra1 and dec1 have different lengths')

    cdef vector[S2Point] vertices1
    cdef int i
    for i in range(n):
        vertices1.push_back(
            S2LatLng.FromRadians(dec1[i], ra1[i])
            .Normalized().ToPoint()
        )

    cdef unique_ptr[S2Loop] loop1
    loop1 = make_unique[S2Loop](vertices1)
    loop1.get().Normalize()
    cdef S2Polygon poly1
    poly1.Init(move(loop1))

    n = ra2.shape[0]
    if dec2.shape[0] != n:
        raise ValueError('ra2 and dec2 have different lengths')

    cdef vector[S2Point] vertices2
    for i in range(n):
        vertices2.push_back(
            S2LatLng.FromRadians(dec2[i], ra2[i])
            .Normalized().ToPoint()
        )

    cdef unique_ptr[S2Loop] loop2
    loop2 = make_unique[S2Loop](vertices2)
    loop2.get().Normalize()
    cdef S2Polygon poly2
    poly2.Init(move(loop2))

    return poly1.Intersects(&poly2)


def polygon_string_intersects_polygon(s, double[:] ra2, double[:] dec2):
    """Test for polygon intersection."""
    coords = np.radians(np.array([c.split(':') for c in s.split(',')], float))
    return polygon_intersects_polygon(coords[:, 0], coords[:, 1], ra2, dec2)


cdef class SpatialIndexer:
    """Spatial indexer.

    Parameters
    ----------
    min_edge_length : double
        Minimum edge length to index, radians.

    """

    cdef S2RegionTermIndexer _indexer


    def __cinit__(self, min_edge_length):
        cdef S2RegionTermIndexer.Options options
        options.set_max_level(kAvgEdge.GetClosestLevel(min_edge_length))
        self._indexer = S2RegionTermIndexer(options)


    def index_points_by_area(self, ra, dec):
        """Index area covered by point(s).
        
        Bounding box is aligned with RA and Dec.


        Parameters
        ----------
        ra, dec : ndarray
            Coordinates defining an area to index (e.g., corners, set of
            stars), radians.


        Returns
        -------
        terms : list of byte strings


        """

        cdef int n = ra.shape[0]
        if dec.shape[0] != n:
            raise ValueError('ra and dec have different lengths')

        cdef S2LatLngRect rect = S2LatLngRect()
        cdef int i
        for i in range(n):
            rect.AddPoint(S2LatLng.FromRadians(dec[i], ra[i]).Normalized())
        
        # use PolarClosure to ensure the entire pole is included, if needed.
        cdef vector[string] terms = self._indexer.GetIndexTerms(
            rect.PolarClosure(), b"")
        return [str(term) for term in terms]

    def index_polygon_string(self, s):
        """Index the polygon described by this string.
        
        
        Parameters
        ----------
        s : str
            Comma-separated RA:Dec pairs in degrees, e.g., "1:1, 1:2, 2:1".
        

        Returns
        -------
        terms : list of byte strings

        """

        coords = np.radians(np.array([c.split(':') for c in s.split(',')], float))
        return self.index_polygon(coords[:, 0], coords[:, 1])

    
    def index_polygon(self, double[:] ra, double[:] dec):
        """Index the given polygon.


        Parameters
        ----------
        ra, dec : ndarray
            Polygon vertices, radians.  Order does not matter.  It is assumed the
            polygon should be smaller than 1/2 the sphere.


        Returns
        -------
        terms : list of byte strings


        """

        cdef int n = ra.shape[0]
        if dec.shape[0] != n:
            raise ValueError('ra and dec have different lengths')

        cdef vector[S2Point] vertices
        cdef int i
        for i in range(n):
            vertices.push_back(S2LatLng.FromRadians(dec[i], ra[i]).Normalized().ToPoint())

        cdef unique_ptr[S2Loop] loop
        loop = make_unique[S2Loop](vertices)
        loop.get().Normalize()
        cdef S2Polygon poly
        poly.Init(move(loop))

        cdef vector[string] terms = self._indexer.GetIndexTerms(poly, b"")
        return [str(term) for term in terms]


    def query_line(self, double[:] ra, double[:] dec):
        """Query terms for a region defined by a multi-segmented line.


        Parameters
        ----------
        ra, dec : ndarray
            Coordinates defining line vertices, radians.


        Returns
        -------
        terms : list of byte strings

        """

        cdef int n = ra.shape[0]
        if dec.shape[0] != n:
            raise ValueError('ra and dec have different lengths')

        cdef vector[S2LatLng] vertices
        cdef int i
        for i in range(n):
            vertices.push_back(S2LatLng.FromRadians(dec[i], ra[i]).Normalized())
        cdef S2Polyline line = S2Polyline(vertices)

        cdef vector[string] terms = self._indexer.GetQueryTerms(line, b"")
        return [str(term) for term in terms]

    def query_polygon(self, double[:] ra, double[:] dec):
        """Query terms for a polygon.


        Parameters
        ----------
        ra, dec : ndarray
            Polygon vertices, radians.  Order does not matter.  It is assumed
            the polygon should be smaller than 1/2 the sphere.


        Returns
        -------
        terms : list of byte strings

        """

        cdef int n = ra.shape[0]
        if dec.shape[0] != n:
            raise ValueError('ra and dec have different lengths')

        cdef vector[S2Point] vertices
        cdef int i
        for i in range(n):
            vertices.push_back(S2LatLng.FromRadians(dec[i], ra[i]).Normalized().ToPoint())

        cdef unique_ptr[S2Loop] loop
        loop = make_unique[S2Loop](vertices)
        loop.get().Normalize()
        cdef S2Polygon poly
        poly.Init(move(loop))

        cdef vector[string] terms = self._indexer.GetQueryTerms(poly, b"")
        return [str(term) for term in terms]

    def query_about_line(self, double[:] ra, double[:] dec, 
                         double[:] a, double[:] b):
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

        cdef Py_ssize_t n = ra.shape[0]
        cdef Py_ssize_t n_dec = dec.shape[0]
        cdef Py_ssize_t n_a = a.shape[0]
        cdef Py_ssize_t n_b = b.shape[0]
        if any([n != n_dec, n_dec != n_a, n_a != n_b]):
            raise ValueError('All arrays must have equal length')

        cdef int i, j
        # for i in range(n):
        #     print(ra[i], dec[i], a[i], b[i])

        # the query region spine is the input line extended by +/-a
        cdef double[:,:] spine = np.empty((n + 2, 2))
        cdef double[:] _ra = np.empty(2 * (n + 2))
        cdef double[:] _dec = np.empty(2 * (n + 2))

        # direction of each vector as position angle, including extension
        cdef double[:] pa = np.empty(n + 2)
        for i in range(n - 1):
            pa[i + 1] = _position_angle(ra[i], dec[i], ra[i + 1], dec[i + 1])
        pa[0] = pa[1]
        pa[n] = pa[n - 1]
        pa[n + 1] = pa[n]

        spine[0, 0], spine[0, 1] = _offset_by(ra[0], dec[0], pa[0], -a[0])
        for i in range(n):
            spine[i + 1, 0] = ra[i]
            spine[i + 1, 1] = dec[i]
        spine[n + 1, 0], spine[n + 1, 1] = _offset_by(
            ra[n - 1], dec[n - 1], pa[n - 1], a[n - 1])

        # extend b to match
        cdef double[:] _b = np.r_[b[0], b, b[-1]]

        # construct query region polygon
        # first half of the region 
        for i in range(n + 2):
            _ra[i], _dec[i] = _offset_by(spine[i, 0], spine[i, 1], pa[i] + M_PI_2, _b[i])

        # second half
        for i in range(n + 2):
            j = n + 1 - i
            _ra[n + 2 + i], _dec[n + 2 + i] = _offset_by(spine[j, 0], spine[j, 1], pa[j] - M_PI_2, _b[j])

        #for i in range(n * 2 + 4):
        #    print(f'{_ra[i]:.3f}, {_dec[i]:.3f}')

        return self.query_polygon(_ra, _dec), np.array(_ra), np.array(_dec)
