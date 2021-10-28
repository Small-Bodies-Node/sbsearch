# Licensed with the 3-clause BSD license.  See LICENSE for details.

import pytest
import numpy as np
from astropy.coordinates import SkyCoord
from ..spatial import (SpatialIndexer,
                       position_angle,
                       offset_by,
                       verify_build_polygon,
                       polygon_intersects_line,
                       polygon_string_intersects_line,
                       PolygonBuildError)


@pytest.mark.parametrize(
    "ra, dec, pa",
    [([0, 0], [0, 1], 0),
     ([0, 1], [0, 0], np.pi / 2),
     ([0, 0], [0, -1], np.pi),
     ([0, -1], [0, 0], -np.pi / 2)]
)
def test_position_angle(ra, dec, pa):
    assert np.isclose(position_angle(ra[0], dec[0], ra[1], dec[1]), pa)


@pytest.mark.parametrize(
    "p1, pa, rho, p2",
    [((0, 0), 0, 1, (0, 1)),
     ((0, 0), np.pi / 2, 1, (1, 0)),
     ((0, 0), np.pi, 1, (0, -1)),
     ((0, 0), 3 / 2 * np.pi, 1, (-1, 0))]
)
def test_offset_by(p1, pa, rho, p2):
    ra, dec = offset_by(p1[0], p1[1], pa, rho)
    assert np.isclose(ra, p2[0])
    assert np.isclose(dec, p2[1])


def test_build_polygon():
    # simple case
    coords = SkyCoord([0, 0, 1, 0], [0, 1, 1, 0], unit='deg')
    assert verify_build_polygon(coords.ra.rad, coords.dec.rad)

    # polygon with a loop
    coords = SkyCoord([0, 10, 10, 5, 5, 4, 4, 11, 11, 0, 0],
                      [0, 0, 10, 5, -5, -5, 5, 12, -1, -1, 0],
                      unit='deg')
    assert verify_build_polygon(coords.ra.rad, coords.dec.rad)

    # bad polygon
    coords = SkyCoord([0, 0, 1], [0, 1, 1], unit='deg')
    with pytest.raises(PolygonBuildError):
        verify_build_polygon(coords.ra.rad, coords.dec.rad)


@pytest.fixture
def indexer():
    return SpatialIndexer(0.01)


# Index and query terms below were generated with:
#   max_cells = 8
#   min_edge_length = 0.01 rad
#   set_min_level(kAvgEdge.GetClosestLevel(0.17))  # 10 deg


EXPECTED_INDEX_TERMS = {
    '10194',
    '1019',
    '101c',
    '101',
    '104',
    '1019c',
    '$101b',
    '101b',
    '101c4',
    '101d',
    '101cc',
    '101ec',
    '101f',
}


class TestSpatialIndexer:
    def test_index_points_by_area(self, indexer):
        coords = SkyCoord([1, 2], [3, 4], unit='deg')
        indices = indexer.index_points_by_area(coords.ra.rad, coords.dec.rad)
        assert set(indices) == EXPECTED_INDEX_TERMS

    def test_index_polygon(self, indexer):
        coords = SkyCoord([1, 2, 2, 1], [3, 3, 4, 4], unit='deg')
        indices = indexer.index_polygon(coords.ra.rad, coords.dec.rad)
        # should be identical to test_index_points_by_area
        assert set(indices) == EXPECTED_INDEX_TERMS

    def test_index_polygon_string(self, indexer):
        indices = indexer.index_polygon_string('1:3, 2:3, 2:4, 1:4')
        # should be identical to test_index_polygon
        assert set(indices) == EXPECTED_INDEX_TERMS

    def test_query_line(self, indexer):
        coords = SkyCoord([1, 2], [3, 4], unit='deg')
        indices = indexer.query_line(coords.ra.rad, coords.dec.rad)
        expected = {
            '$101',
            '$1019',
            '$101b',
            '$101c',
            '$101d',
            '$104',
            '10194',
            '1019c',
            '101bc',
            '101c4',
            '101cc',
        }
        assert set(indices) == expected

        # should have some values in common with test_index_polygon
        assert len(set(indices) & EXPECTED_INDEX_TERMS) > 0

    def test_query_polygon(self, indexer):
        coords = SkyCoord([1, 2, 2, 1], [3, 3, 4, 4], unit='deg')
        indices = indexer.query_polygon(coords.ra.rad, coords.dec.rad)
        expected = {
            '$101',
            '$1019',
            '$101c',
            '$101d',
            '$101f',
            '$104',
            '10194',
            '1019c',
            '101b',
            '101c4',
            '101cc',
            '101ec',
        }
        assert set(indices) == expected

        # reverse the polygon, should be the same result
        coords = SkyCoord([1, 1, 2, 2], [3, 4, 4, 3], unit='deg')
        indices = indexer.query_polygon(coords.ra.rad, coords.dec.rad)
        assert set(indices) == expected

    @pytest.mark.parametrize(
        "ra, dec, a, b, n",
        [
            ((0.0, 0.0), (0.0, 0.01), (0.01, 0.01), (0.01, 0.01), 0),
            ((0.01745329, 0.03490659), (0.05235988, 0.06981317),
             (0.01, 0.01), (0.01, 0.01), 4),
            ((0.0, 0.01745329), (0.0, 0.05235988), (0.01, 0.01), (0.01, 0.01), 2),
        ]
    )
    def test_query_about_line(self, ra, dec, a, b, n, indexer):
        ra = np.array(ra)
        dec = np.array(dec)
        a = np.array(a)
        b = np.array(b)
        indices, poly_ra, poly_dec = indexer.query_about_line(ra, dec, a, b)

        assert len(poly_ra) == len(poly_dec)
        assert polygon_intersects_line(poly_ra, poly_dec, ra, dec)

        # should have n values in common with test_index_polygon
        assert len(set(indices) & EXPECTED_INDEX_TERMS) == n

    @pytest.mark.parametrize(
        "poly_ra, poly_dec, line_ra, line_dec, intersects",
        [
            ((0, 1, 1), (0, 0, 1), (0, 1), (0, 0.25), True),
            ((0, 1, 1), (0, 0, 1), (0, -1), (0, 0.25), False),
        ]
    )
    def test_polygon_intersects_line(self, poly_ra, poly_dec, line_ra,
                                     line_dec, intersects):
        assert polygon_intersects_line(
            np.array(poly_ra, float), np.array(poly_dec, float),
            np.array(line_ra, float), np.array(line_dec, float)
        ) == intersects

    @pytest.mark.parametrize(
        "s, line_ra, line_dec, intersects",
        [
            ('0:0, 1:0, 1:1', (0, 1), (0, 0.25), True),
            ('0:0, 1:0, 1:1', (0, -1), (0, 0.25), False),
        ]
    )
    def test_polygon_string_intersects_line(self, s, line_ra,
                                            line_dec, intersects):
        assert polygon_string_intersects_line(
            s, np.array(line_ra, float), np.array(line_dec, float)
        ) == intersects

# def test_S2PointOrientation():
#     '''Verify that x,y,z==0,0,1==ra,dec==*,+90 deg'''
#     assert np.isclose(_test_turn_angle(), np.pi / 2)
# class TestAngle:
#     def test_init(self):
#         a: Angle = Angle()
#         assert a.radians() == 0
#     def test_init_radians(self):
#         a: Angle = Angle.Radians(1)
#         assert a.radians() == 1
#     def test_init_degrees(self):
#         a: Angle = Angle.Degrees(1)
#         assert a.degrees() == 1
#     def test_degrees(self):
#         a: Angle = Angle.Radians(1)
#         assert np.isclose(a.degrees(), np.degrees(1))
#     def test_radians(self):
#         a: Angle = Angle.Degrees(1)
#         assert np.isclose(a.radians(), np.radians(1))
#     def test_Normalized(self):
#         a: Angle = Angle.Degrees(200)
#         b: Angle = a.Normalized()
#         assert np.isclose(a.degrees(), 200)
#         assert np.isclose(b.degrees(), -160)
#     def test_Normalize(self):
#         a: Angle = Angle.Degrees(200)
#         a.Normalize()
#         assert np.isclose(a.degrees(), -160)
# class TestPoint:
#     def test_init(self):
#         p: Point = Point(0, 1, 2)
#         assert p.x == 0
#         assert p.y == 1
#         assert p.z == 2
#     def test_getitem(self):
#         p: Point = Point(0, 1, 2)
#         assert p[0] == 0
#         assert p[1] == 1
#         assert p[2] == 2
