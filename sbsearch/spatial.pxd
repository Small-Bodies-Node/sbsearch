# Licensed with the 3-clause BSD license.  See LICENSE for details.
# cython: language_level=3

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr
from libcpp cimport bool

cdef extern from "libspatial.cpp":
    double _position_angle(double, double, double, double)
    void _offset_by(double, double, double, double, double&, double&)
    S2Point _offset_point_by(S2Point point, double pa, double rho)
    void _build_polygon(double*, double*, int, bool, S2Polygon&)
    void _build_polygon_from_vertices(vector[S2Point], bool, S2Polygon&)
    S2Polyline _build_polyline(double*, double*, int, double, double)
    bool _verify_build_polygon(double*, double*, int) except +
    bool _polygon_contains_point(double*, double*, int, double, double)
    bool _polygon_intersects_line(double*, double*, int, double*, double*, int, double, double)
    bool _polygon_intersects_about_line(double*, double*, int, double*, double*, int, double*, double*, double, double)
    bool _polygon_intersects_polygon(double*, double*, int, double*, double*, int)
    bool _polygon_intersects_cap(double*, double*, int, double, double, double, int)

cdef extern from "s2/s2metrics.h" namespace "S2":
    cdef cppclass LengthMetric:
        int GetClosestLevel(double)

    cdef LengthMetric kAvgEdge

cdef extern from "s2/s2point.h":
    cdef cppclass S2Point:
        S2Point()
        S2Point(double, double, double)
        float x()
        float y()
        float z()

cdef extern from "s2/s2shape.h":
    cdef cppclass S2Shape:
        S2Shape()

cdef extern from "s2/s2region.h":
    cdef cppclass S2Region:
        S2Region()

cdef extern from "s2/s2polygon.h":
    cdef cppclass S2Polygon(S2Region):
        S2Polygon()
        S2Polygon(unique_ptr[S2Loop])
        void Init(unique_ptr[S2Loop])
        bool IsValid()
        bool Intersects(const S2Polyline&)
        bool Intersects(const S2Polygon*)
        bool Contains(const S2Point&)

        cppclass Shape(S2Shape):
            Shape()

cdef extern from "s2/s2loop.h":
    cdef cppclass S2Loop(S2Region):
        S2Loop()
        S2Loop(const vector[S2Point]&)
        void Normalize()

cdef extern from "s2/s2polyline.h":
    cdef cppclass S2Polyline(S2Region):
        S2Polyline()
        S2Polyline(const vector[S2LatLng]&)
        S2Polyline(const vector[S2Point]&)
        S2Point Interpolate(double)
        S2Point GetSuffix(double, int*)
        S2Point& vertex(int)

cdef extern from "s2/s2latlng.h":
    cdef cppclass S2LatLng:
        S2LatLng()
        S2LatLng(S2LatLng&)
        S2LatLng(const S2Point&)
        S2LatLng(S1Angle, S1Angle)
        S2Point ToPoint()
        @staticmethod
        S2LatLng FromRadians(double, double)
        S2LatLng Normalized()
        S1Angle lat()
        S1Angle lng()

cdef extern from "s2/s1angle.h":
    cdef cppclass S1Angle:
        S1Angle()
        @staticmethod
        S1Angle Radians(double)
        @staticmethod
        S1Angle Degrees(double)
        double radians()
        double degrees()
