# Licensed with the 3-clause BSD license.  See LICENSE for details.
# cython: language_level=3

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr
from libcpp cimport bool

cdef extern from "s2/s2metrics.h" namespace "S2":
    cdef cppclass LengthMetric:
        int GetClosestLevel(double)

    cdef LengthMetric kAvgEdge

cdef extern from "s2/s1angle.h":
    cdef cppclass S1Angle:
        S1Angle()
        @staticmethod
        S1Angle Radians(double)
        @staticmethod
        S1Angle Degrees(double)
        double radians()
        double degrees()

cdef extern from "s2/s2point.h":
    cdef cppclass S2Point:
        S2Point()
        S2Point(double, double, double)

cdef extern from "s2/s2latlng.h":
    cdef cppclass S2LatLng:
        S2LatLng()
        S2LatLng(S2LatLng&)
        S2LatLng(S1Angle, S1Angle)
        S2Point ToPoint()
        @staticmethod
        S2LatLng FromRadians(double, double)
        S2LatLng Normalized()
        
cdef extern from "s2/s2region.h":
    cdef cppclass S2Region:
        S2Region()

cdef extern from "s2/s2latlng_rect.h":
    cdef cppclass S2LatLngRect(S2Region):
        S2LatLngRect()
        S2LatLngRect(const S2LatLng&, const S2LatLng&)
        void AddPoint(const S2LatLng&)
        S2LatLngRect PolarClosure()

cdef extern from "s2/s2latlng_rect_bounder.h":
    cdef cppclass S2LatLngRectBounder:
        S2LatLngRectBounder()
        void AddLatLng(const S2LatLng&)
        S2LatLngRect GetBound()

cdef extern from "s2/s2polyline.h":
    cdef cppclass S2Polyline(S2Region):
        S2Polyline()
        S2Polyline(const vector[S2LatLng]&)
        S2Polyline(const vector[S2Point]&)
        S2Point Interpolate(double)
        S2Point GetSuffix(double, int*)

cdef extern from "s2/s2region_term_indexer.h":
    cdef cppclass S2RegionTermIndexer:
        cppclass Options:
            void set_max_level(int)
            void set_min_level(int)

        S2RegionTermIndexer()
        S2RegionTermIndexer(const Options&)

        vector[string] GetIndexTerms(const S2Region&, string)
        vector[string] GetQueryTerms(const S2Region&, string)

cdef extern from "s2/s2loop.h":
    cdef cppclass S2Loop(S2Region):
        S2Loop()
        S2Loop(const vector[S2Point]&)
        void Normalize()

cdef extern from "s2/s2debug.h":
    cdef cppclass S2Debug:
        pass

    cdef S2Debug ALLOW
    cdef S2Debug DISABLE

cdef extern from "s2/s2polygon.h":
    cdef cppclass S2Polygon(S2Region):
        S2Polygon()
        S2Polygon(unique_ptr[S2Loop])
        void Init(unique_ptr[S2Loop])
        bool Intersects(const S2Polyline&)
        bool Intersects(const S2Polygon*)