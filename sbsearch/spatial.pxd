# Licensed with the 3-clause BSD license.  See LICENSE for details.
# cython: language_level=3

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr
from libcpp cimport bool

cdef extern from "libspatial.cpp":
    enum IntersectionType:
        PolygonContainsPoint = 0,
        PolygonContainsArea = 1,
        PolygonIntersectsArea = 2,
        AreaContainsPolygon = 3

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
    bool _polygon_intersects_cap(double*, double*, int, double, double, double, IntersectionType)

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

cdef extern from "s2/s1chord_angle.h":
    cdef cppclass S1ChordAngle:
        @staticmethod
        S1ChordAngle Radians(double)

        bool operator>(S1ChordAngle, S1ChordAngle)
        bool operator>=(S1ChordAngle, S1ChordAngle)
        bool operator<(S1ChordAngle, S1ChordAngle)
        bool operator<=(S1ChordAngle, S1ChordAngle)

cdef extern from "s2/s2point.h":
    cdef cppclass S2Point:
        S2Point()
        S2Point(double, double, double)
        float x()
        float y()
        float z()

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
        S2Point& vertex(int)

cdef extern from "s2/s2region_term_indexer.h":
    cdef cppclass S2RegionTermIndexer:
        cppclass Options:
            void set_max_level(int)
            void set_min_level(int)
            int max_cells()
            void set_max_cells(int)

        S2RegionTermIndexer()
        S2RegionTermIndexer(const Options&)

        vector[string] GetIndexTerms(const S2Region&, absl_string_view)
        vector[string] GetQueryTerms(const S2Region&, absl_string_view)
        vector[string] GetQueryTerms(const S2Point&, absl_string_view)

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

cdef extern from "s2/s2shape.h":
    cdef cppclass S2Shape:
        S2Shape()

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

cdef extern from "s2/s2error.h" namespace "S2Error":
    cdef enum Code:
        OK = 0,  # // No error.

        # ////////////////////////////////////////////////////////////////////
        # // Generic errors, not specific to geometric objects:

        UNKNOWN = 1000,              # // Unknown error.
        UNIMPLEMENTED = 1001,        # // Operation is not implemented.
        OUT_OF_RANGE = 1002,         # // Argument is out of range.
        INVALID_ARGUMENT = 1003,     # // Invalid argument (other than a range error).
        FAILED_PRECONDITION = 1004,  # // Object is not in the required state.
        INTERNAL = 1005,             # // An internal invariant has failed.
        DATA_LOSS = 1006,            # // Data loss or corruption.
        RESOURCE_EXHAUSTED = 1007,   # // A resource has been exhausted.

        # ////////////////////////////////////////////////////////////////////
        # // Error codes in the following range can be defined by clients:

        USER_DEFINED_START = 1000000,
        USER_DEFINED_END   = 9999999,

        # ////////////////////////////////////////////////////////////////////
        # // Errors that apply to more than one type of geometry:

        NOT_UNIT_LENGTH = 1,     # // Vertex is not unit length.
        DUPLICATE_VERTICES = 2,  # // There are two identical vertices.
        ANTIPODAL_VERTICES = 3,  # // There are two antipodal vertices.

        # ////////////////////////////////////////////////////////////////////
        # // S2Loop errors:

        LOOP_NOT_ENOUGH_VERTICES = 100,  # // Loop with fewer than 3 vertices.
        LOOP_SELF_INTERSECTION = 101,    # // Loop has a self-intersection.

        # ////////////////////////////////////////////////////////////////////
        # // S2Polygon errors:

        POLYGON_LOOPS_SHARE_EDGE = 200,  # // Two polygon loops share an edge.
        POLYGON_LOOPS_CROSS = 201,       # // Two polygon loops cross.
        POLYGON_EMPTY_LOOP = 202,        # // Polygon has an empty loop.
        POLYGON_EXCESS_FULL_LOOP = 203,  # // Non-full polygon has a full loop.

        # // InitOriented() was called and detected inconsistent loop orientations.
        POLYGON_INCONSISTENT_LOOP_ORIENTATIONS = 204,

        # // Loop depths don't correspond to any valid nesting hierarchy.
        POLYGON_INVALID_LOOP_DEPTH = 205,

        # // Actual polygon nesting does not correspond to the nesting hierarchy
        # // encoded by the loop depths.
        POLYGON_INVALID_LOOP_NESTING = 206,

        # ////////////////////////////////////////////////////////////////////
        # // S2Builder errors:

        # // The S2Builder snap function moved a vertex by more than the specified
        # // snap radius.
        BUILDER_SNAP_RADIUS_TOO_SMALL = 300,

        # // S2Builder expected all edges to have siblings (as specified by
        # // S2Builder::GraphOptions::SiblingPairs::REQUIRE), but some were missing.
        BUILDER_MISSING_EXPECTED_SIBLING_EDGES = 301,

        # // S2Builder found an unexpected degenerate edge.  For example,
        # // Graph::GetLeftTurnMap() does not support degenerate edges.
        BUILDER_UNEXPECTED_DEGENERATE_EDGE = 302,

        # // S2Builder found a vertex with (indegree != outdegree), which means
        # // that the given edges cannot be assembled into loops.
        BUILDER_EDGES_DO_NOT_FORM_LOOPS = 303,

        # // The edges provided to S2Builder cannot be assembled into a polyline.
        BUILDER_EDGES_DO_NOT_FORM_POLYLINE = 304,

        # // There was an attempt to assemble a polygon from degenerate geometry
        # // without having specified a predicate to decide whether the output is
        # // the empty polygon (containing no points) or the full polygon
        # // (containing all points).
        BUILDER_IS_FULL_PREDICATE_NOT_SPECIFIED = 305
