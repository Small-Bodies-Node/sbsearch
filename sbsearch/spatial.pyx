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
    IntersectionType,
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
    return kAvgEdge.GetClosestLevel(edge_length)


def position_angle(double ra1, double dec1, double ra2, double dec2):
    return _position_angle(ra1, dec1, ra2, dec2)


def offset_by(double ra, double dec, double pa, double rho):
    cdef double new_ra = 0, new_dec = 0
    _offset_by(ra, dec, pa, rho, new_ra, new_dec)
    return new_ra, new_dec


# cdef build_polygon(double[:] ra, double[:] dec, S2Polygon& polygon,
#                     close=True):
#     """Build polygon from vertices.

#     Uses S2Builder in order to accomodate loops.

#     The vertices must form a closed shape, e.g., last vertex = first vertex.
#     Use close=True if they do not.

#     """

#     cdef int n = len(ra)
  
#     cdef vector[S2Point] vertices
#     cdef int i
#     for i in range(n):
#         vertices.push_back(
#             S2LatLng.FromRadians(dec[i], ra[i])
#             .Normalized().ToPoint()
#         )

#     cdef S2Builder.Options builder_options
#     builder_options.set_split_crossing_edges(True)
#     cdef S2Builder* builder = new S2Builder(builder_options)

#     # cdef S2Polygon polygon
#     cdef S2PolygonLayer.Options layer_options
#     layer_options.set_edge_type(EdgeType.UNDIRECTED)
#     builder.StartLayer(make_unique[S2PolygonLayer](&polygon, layer_options))
#     for i in range(1, n):
#         builder.AddEdge(vertices[i - 1], vertices[i])
#     if close:
#         builder.AddEdge(vertices[n - 1], vertices[0])

#     cdef S2Error error
#     builder.Build(&error)
#     if not error.ok():
#         raise PolygonBuildError('({}): {}'.format(
#             error.code(), error.text().decode()))

# def verify_build_polygon(double[:] ra, double[:] dec):
#     """Python interface with build_polygon for testing purposes.
    
#     Does not close the polygon.
    
#     """

#     cdef S2Polygon polygon
#     build_polygon(ra, dec, polygon, close=False)
#     return polygon.IsValid()

# def polygon_contains_point(double[:] poly_ra, double[:] poly_dec,
#                            double point_ra, double point_dec):
#     """Test if the polygon covers the point.
    

#     Parameters
#     ----------
#     poly_ra, poly_dec : ndarray
#         Polygon RA and Dec, radians.

#     point_ra, point_dec : ndarray
#         Point RA and Dec, radians.


#     Returns
#     -------
#     intersects : bool
    
#     """

#     # cdef S2Polygon polygon
#     polygon = s2.S2Polygon()
#     build_polygon(poly_ra, poly_dec, polygon)

#     point = (s2.S2LatLng.FromRadians(point_dec, point_ra)
#                         .Normalized().ToPoint())

#     return polygon.Contains(point)


# def polygon_intersects_line(double[:] poly_ra, double[:] poly_dec,
#                             double[:] line_ra, double[:] line_dec,
#                             double line_start=0, double line_stop=1):
#     """Test if the polygon intersects the line.
    

#     Parameters
#     ----------
#     poly_ra, poly_dec : ndarray
#         Polygon RA and Dec, radians.

#     line_ra, line_dec : ndarray
#         Line RA and Dec, radians.

#     line_start, line_stop : float
#         Test a sub-portion of the line, indicated by the given line fractions
#         (linear interpolation).


#     Returns
#     -------
#     intersects : bool
    
#     """

#     cdef int n = poly_ra.shape[0]
#     if poly_dec.shape[0] != n:
#         raise ValueError('ra and dec have different lengths')

    # cdef S2Polygon polygon
    # build_polygon(poly_ra, poly_dec, polygon)

    # n = line_ra.shape[0]
    # if line_dec.shape[0] != n:
    #     raise ValueError('ra and dec have different lengths')

    # cdef vector[S2Point] line_vertices
    # for i in range(n):
    #     line_vertices.push_back(
    #         S2LatLng.FromRadians(line_dec[i], line_ra[i])
    #         .Normalized().ToPoint()
    #     )
    # cdef S2Polyline line = S2Polyline(line_vertices)

    # # interpolate to a sub-set?
    # cdef S2Point start = line.Interpolate(line_start)
    # cdef S2Point stop = line.Interpolate(line_stop)
    # cdef int next_vertex, last_vertex
    # cdef vector[S2Point] new_line_vertices
    # if line_start != 0 or line_stop != 1:
    #     if line_stop <= line_start:
    #         raise ValueError('line_stop <= line_start')

    #     line.GetSuffix(line_start, &next_vertex)
    #     line.GetSuffix(line_stop, &last_vertex)
    #     # build the new line
    #     new_line_vertices.push_back(start)
    #     if next_vertex != last_vertex:
    #         for i in range(last_vertex, next_vertex):
    #             new_line_vertices.push_back(line_vertices[i])
    #     new_line_vertices.push_back(stop)

    #     line = S2Polyline(new_line_vertices)
   
    # return polygon.Intersects(line)

def verify_build_polygon(ra, dec):
    cdef double[:] _ra = ra
    cdef double[:] _dec = dec
    try:
        return _verify_build_polygon(&_ra[0], &_dec[0], len(ra))
    except ValueError as exc:
        raise PolygonBuildError(str(exc))


def polygon_intersects_line(double[:] poly_ra, double[:] poly_dec,
                            double[:] line_ra, double[:] line_dec,
                            double line_start=0, double line_stop=1):
    return _polygon_intersects_line(&poly_ra[0], &poly_dec[0], len(poly_ra),
                                    &line_ra[0], &line_dec[0], len(line_ra),
                                    line_start, line_stop)

def polygon_string_intersects_line(s, double[:] line_ra, double[:] line_dec,
                                   line_start=0, line_stop=1):
    """Test if the polygon intersects line.
    

    Parameters
    ----------
    s : str
        Comma-separated RA:Dec pairs in degrees, e.g., "1:1, 1:2, 2:1".
    
    line_ra, line_dec : ndarray
        Line RA and Dec, radians.

    line_start, line_stop : float
        Test a sub-portion of the line, indicated by the given line fractions
        (linear interpolation).

    """

    if len(line_ra) != len(line_dec):
        raise ValueError

    coords = np.radians(np.array([c.split(':') for c in s.split(',')], float)).T
    cdef double[:] poly_ra = coords[0]
    cdef double[:] poly_dec = coords[1]
    return _polygon_intersects_line(&poly_ra[0], &poly_dec[0], len(coords),
                                   &line_ra[0], &line_dec[0], len(line_ra),
                                   line_start, line_stop)


# def polygon_intersects_about_line(double[:] poly_ra, double[:] poly_dec,
#                                   double[:] line_ra, double[:] line_dec,
#                                   double[:] a, double[:] b,
#                                   double line_start=0, double line_stop=1):
#     """Test for intersection between polygon and region about a line.
    

#     Parameters
#     ----------
#     poly_ra, poly_dec : ndarray
#         Polygon RA and Dec, radians.

#     line_ra, line_dec : ndarray
#         Line RA and Dec, radians.

#     line_start, line_stop : float
#         Test a sub-portion of the line, indicated by the given line fractions.

#     a, b : ndarray
#         Angular size of the region: ``a`` is along the line, ``b`` 
#         is perpendicular to it, radians.  Only the first and last
#         ``a`` are considered.

#     line_start, line_stop : float
#         Test a sub-portion of the line, indicated by the given line fractions
#         (linear interpolation).


#     Returns
#     -------
#     intersects : bool

#     """

#     cdef int i, j

#     # validate line inputs; polygon is validated in polygon_intersects_polygon
#     cdef int n = line_ra.shape[0]
#     if any([n != line_dec.shape[0], n != a.shape[0], n != b.shape[0]]):
#         raise ValueError('Line arrays must all have equal length')

#     # define the line or sub-segment thereof
#     cdef vector[S2Point] line_vertices
#     for i in range(n):
#         line_vertices.push_back(
#             S2LatLng.FromRadians(line_dec[i], line_ra[i])
#             .Normalized().ToPoint()
#         )
#     cdef S2Polyline line = S2Polyline(line_vertices)

#     # interpolate to a sub-segment?
#     cdef S2Point start = line.Interpolate(line_start)
#     cdef S2Point stop = line.Interpolate(line_stop)
#     cdef int next_vertex, last_vertex
#     cdef vector[S2Point] new_line_vertices
#     if line_start != 0 or line_stop != 1:
#         if line_stop <= line_start:
#             raise ValueError('line_stop <= line_start')

#         line.GetSuffix(line_start, &next_vertex)
#         line.GetSuffix(line_stop, &last_vertex)
#         # build the new line
#         new_line_vertices.push_back(start)
#         if next_vertex != last_vertex:
#             for i in range(last_vertex, next_vertex):
#                 new_line_vertices.push_back(line_vertices[i])
#         new_line_vertices.push_back(stop)

#         line = S2Polyline(new_line_vertices)

#     # transform the line into a polygon

#     # the spine is the input line extended by +/-a
#     # the polygon is the region betweeh spine+b and spine-b
#     cdef double[:,:] spine = np.empty((n + 2, 2))
#     cdef double[:] _ra = np.empty(2 * (n + 2))
#     cdef double[:] _dec = np.empty(2 * (n + 2))

#     # direction of each vector as position angle, including extension
#     cdef double[:] pa = np.empty(n + 2)
#     for i in range(n - 1):
#         pa[i + 1] = _position_angle(line_ra[i], line_dec[i], line_ra[i + 1], line_dec[i + 1])
#     pa[0] = pa[1]
#     pa[n] = pa[n - 1]
#     pa[n + 1] = pa[n]

#     spine[0, 0], spine[0, 1] = _offset_by(line_ra[0], line_dec[0], pa[0], -a[0])
#     for i in range(n):
#         spine[i + 1, 0] = line_ra[i]
#         spine[i + 1, 1] = line_dec[i]
#     spine[n + 1, 0], spine[n + 1, 1] = _offset_by(
#         line_ra[n - 1], line_dec[n - 1], pa[n], a[n - 1])

#     # pad out b to match the number of spine vertices
#     cdef double[:] _b = np.r_[b[0], b, b[n - 1]]

#     # construct polygon
#     # first half of the region 
#     for i in range(n + 2):
#         _ra[i], _dec[i] = _offset_by(spine[i, 0], spine[i, 1], pa[i] + M_PI_2, _b[i])

#     # second half
#     for i in range(n + 2):
#         j = n + 1 - i
#         _ra[n + 2 + i], _dec[n + 2 + i] = _offset_by(spine[j, 0], spine[j, 1], pa[j] - M_PI_2, _b[j])

#     return polygon_intersects_polygon(poly_ra, poly_dec, _ra, _dec)


# def polygon_intersects_polygon(double[:] ra1, double[:] dec1,
#                                double[:] ra2, double[:] dec2):
#     """Test for polygon intersection."""
#     cdef int n = ra1.shape[0]
#     if dec1.shape[0] != n:
#         raise ValueError('ra1 and dec1 have different lengths')

#     n = ra2.shape[0]
#     if dec2.shape[0] != n:
#         raise ValueError('ra2 and dec2 have different lengths')

#     cdef S2Polygon polygon1, polygon2
#     build_polygon(ra1, dec1, polygon1)
#     build_polygon(ra2, dec2, polygon2)

#     return polygon1.Intersects(&polygon2)

def polygon_intersects_about_line(double[:] poly_ra, double[:] poly_dec,
                                  double[:] line_ra, double[:] line_dec,
                                  double[:] a, double[:] b,
                                  double line_start=0, double line_stop=1):
    return _polygon_intersects_about_line(&poly_ra[0], &poly_dec[0], len(poly_ra),
                                          &line_ra[0], &line_dec[0], len(line_ra),
                                          &a[0], &b[0], line_start, line_stop)


def polygon_string_intersects_about_line(s, double[:] line_ra, double[:] line_dec,
                                         double[:] a, double[:] b, line_start=0, line_stop=1):
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

    coords = np.radians(np.array([c.split(':') for c in s.split(',')], float)).T
    cdef double[:] poly_ra = coords[0]
    cdef double[:] poly_dec = coords[1]
    return polygon_intersects_about_line(poly_ra, poly_dec, line_ra, line_dec,
                                         a, b, line_start, line_stop)


cdef polygon_intersects_polygon(double[:] ra1, double[:] dec1, double[:] ra2, double[:] dec2):
    return _polygon_intersects_polygon(&ra1[0], &dec1[0], len(ra1), &ra2[0], &dec2[0], len(ra2))


def polygon_string_intersects_polygon(s, double[:] ra2, double[:] dec2):
    """Test for polygon intersection."""
    coords = np.radians(np.array([c.split(':') for c in s.split(',')], float)).T
    return polygon_intersects_polygon(coords[0], coords[1], ra2, dec2)


def polygon_string_intersects_polygon_string(s1, s2):
    """Test for polygon intersection."""
    coords1 = np.radians(np.array([c.split(':') for c in s1.split(',')], float))
    coords2 = np.radians(np.array([c.split(':') for c in s2.split(',')], float))
    return polygon_intersects_polygon(coords1[0], coords1[1], coords2[0], coords2[1])


def polygon_contains_point(double[:] poly_ra, double[:] poly_dec, double point_ra, double point_dec):
    return _polygon_contains_point(&poly_ra[0], &poly_dec[0], len(poly_ra), point_ra, point_dec)


def polygon_string_contains_point(s, double point_ra, double point_dec):
    """Test for polygon intersection."""
    coords = np.radians(np.array([c.split(':') for c in s.split(',')], float)).T
    return polygon_contains_point(coords[0], coords[1], point_ra, point_dec)


def polygon_intersects_cap(double[:] poly_ra, double[:] poly_dec, double point_ra, double point_dec,
                           double radius, int intersection_type):
    return _polygon_intersects_cap(&poly_ra[0], &poly_dec[0], len(poly_ra), point_ra, point_dec,
                                   radius, IntersectionType(intersection_type))


def polygon_string_intersects_cap(s, double point_ra, double point_dec, double radius, int intersection_type):
    coords = np.radians(np.array([c.split(':') for c in s.split(',')], float)).T
    return polygon_intersects_cap(coords[0], coords[1], point_ra, point_dec, radius, intersection_type)


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

    latlng = s2.S2LatLng()

    ra, dec = np.empty((2, 4))
    for i in range(4):
        latlng = s2.S2LatLng(cell.GetVertex(i))
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

    def index_polygon_string(self, s):
        """Index the polygon described by this string.
        
        
        Parameters
        ----------
        s : str
            Comma-separated RA:Dec pairs in degrees, e.g., "1:1, 1:2, 2:1".
        

        Returns
        -------
        terms : list of strings

        """

        coords = np.radians(np.array([c.split(':') for c in s.split(',')], float))
        return self.index_polygon(coords[:, 0], coords[:, 1])

    @staticmethod
    def vertices_to_loop(ra, dec):
        loop = s2.S2Loop()
        loop.Init([s2.S2LatLng.FromRadians(*v).ToPoint() for v in zip(dec, ra)])
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
            The point to query.


        Returns
        -------
        terms : list of strings
        
        """
        
        p = s2.S2LatLng.FromRadians(dec, ra).ToPoint()
        terms = self._indexer.GetQueryTerms(p, "")
        return terms

    @staticmethod
    def vertices_to_polyline(ra, dec):
        vertices = [s2.S2LatLng.FromRadians(*v) for v in zip(dec, ra)]
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