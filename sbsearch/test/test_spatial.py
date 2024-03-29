# Licensed with the 3-clause BSD license.  See LICENSE for details.

import pytest
import numpy as np
from astropy.coordinates import SkyCoord
from ..core import polygon_string_to_arrays
from ..spatial import (
    SpatialIndexer,
    position_angle,
    offset_by,
    verify_build_polygon,
    polygon_intersects_line,
    polygon_intersects_line,
    polygon_intersects_polygon,
    polygon_intersects_cap,
    PolygonBuildError,
)


@pytest.mark.parametrize(
    "ra, dec, pa",
    [
        ([0, 0], [0, 1], 0),
        ([0, 1], [0, 0], np.pi / 2),
        ([0, 0], [0, -1], np.pi),
        ([0, -1], [0, 0], -np.pi / 2),
    ],
)
def test_position_angle(ra, dec, pa):
    assert np.isclose(position_angle(ra[0], dec[0], ra[1], dec[1]), pa)


@pytest.mark.parametrize(
    "p1, pa, rho, p2",
    [
        ((0, 0), 0, 1, (0, 1)),
        ((0, 0), np.pi / 2, 1, (1, 0)),
        ((0, 0), np.pi, 1, (0, -1)),
        ((0, 0), 3 / 2 * np.pi, 1, (-1, 0)),
    ],
)
def test_offset_by(p1, pa, rho, p2):
    ra, dec = offset_by(p1[0], p1[1], pa, rho)
    assert np.isclose(ra, p2[0])
    assert np.isclose(dec, p2[1])


def test_build_polygon():
    # simple case
    coords = SkyCoord([0, 0, 1, 0], [0, 1, 1, 0], unit="deg")
    assert verify_build_polygon(coords.ra.rad, coords.dec.rad)

    # polygon with a loop
    coords = SkyCoord(
        [0, 10, 10, 5, 5, 4, 4, 11, 11, 0, 0],
        [0, 0, 10, 5, -5, -5, 5, 12, -1, -1, 0],
        unit="deg",
    )
    assert verify_build_polygon(coords.ra.rad, coords.dec.rad)

    # bad polygon
    coords = SkyCoord([0, 0, 1], [0, 1, 1], unit="deg")
    with pytest.raises(PolygonBuildError):
        verify_build_polygon(coords.ra.rad, coords.dec.rad)


@pytest.fixture
def indexer():
    # Index and query terms below were generated with:
    #   max_cells = 8
    #   min_edge_length = 0.01 rad
    #   max_edge_length = 0.17 rad
    return SpatialIndexer(0.01, 0.17)


EXPECTED_INDEX_TERMS = {
    "10194",
    "1019",
    "101c",
    "101",
    "104",
    "1019c",
    "$101b",
    "101b",
    "101c4",
    "101d",
    "101cc",
    "101ec",
    "101f",
}


class TestSpatialIndexer:
    def test_index_polygon(self, indexer):
        coords = SkyCoord([1, 2, 2, 1], [3, 3, 4, 4], unit="deg")
        indices = indexer.index_polygon(coords.ra.rad, coords.dec.rad)
        assert set(indices) == EXPECTED_INDEX_TERMS

    def test_query_cap(self, indexer):
        expected = {  # similar to query polygon below, but a few more terms
            "$101",
            "101c4",
            "101b4",
            "101a4",
            "101bc",
            "$104",
            "$101f",
            "10194",
            "101cc",
            "$101b",
            "$1019",
            "$101d",
            "$101c",
            "1019c",
            "101ec",
        }
        indices = indexer.query_cap(np.radians(1.5), np.radians(3.5), np.radians(0.5))
        assert set(indices) == expected

        # should have some values in common with test_index_polygon
        assert len(set(indices) & EXPECTED_INDEX_TERMS) > 0

    def test_query_line(self, indexer):
        coords = SkyCoord([1, 2], [3, 4], unit="deg")
        indices = indexer.query_line(coords.ra.rad, coords.dec.rad)
        expected = {
            "$101",
            "$1019",
            "$101b",
            "$101c",
            "$101d",
            "$104",
            "10194",
            "1019c",
            "101bc",
            "101c4",
            "101cc",
        }
        assert set(indices) == expected

        # should have some values in common with test_index_polygon
        assert len(set(indices) & EXPECTED_INDEX_TERMS) > 0

    def test_query_polygon(self, indexer):
        coords = SkyCoord([1, 2, 2, 1], [3, 3, 4, 4], unit="deg")
        indices = indexer.query_polygon(coords.ra.rad, coords.dec.rad)
        expected = {
            "$101",
            "$1019",
            "$101c",
            "$101d",
            "$101f",
            "$104",
            "10194",
            "1019c",
            "101b",
            "101c4",
            "101cc",
            "101ec",
        }
        assert set(indices) == expected

        # reverse the polygon, should be the same result
        coords = SkyCoord([1, 1, 2, 2], [3, 4, 4, 3], unit="deg")
        indices = indexer.query_polygon(coords.ra.rad, coords.dec.rad)
        assert set(indices) == expected

    @pytest.mark.parametrize(
        "ra, dec, a, b, n",
        [
            ((0.0, 0.0), (0.0, 0.01), (0.01, 0.01), (0.01, 0.01), 0),
            (
                (0.01745329, 0.03490659),
                (0.05235988, 0.06981317),
                (0.01, 0.01),
                (0.01, 0.01),
                4,
            ),
            ((0.0, 0.01745329), (0.0, 0.05235988), (0.01, 0.01), (0.01, 0.01), 2),
        ],
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
    ],
)
def test_polygon_intersects_line(poly_ra, poly_dec, line_ra, line_dec, intersects):
    assert (
        polygon_intersects_line(
            np.radians(poly_ra),
            np.radians(poly_dec),
            np.radians(line_ra),
            np.radians(line_dec),
        )
        == intersects
    )


@pytest.mark.parametrize(
    "s, line_ra, line_dec, intersects",
    [
        ("0:0, 1:0, 1:1", (0, 1), (0, 0.25), True),
        ("0:0, 1:0, 1:1", (0, -1), (0, 0.25), False),
    ],
)
def test_polygon_intersects_line(s, line_ra, line_dec, intersects):
    poly_ra, poly_dec = polygon_string_to_arrays(s)
    assert (
        polygon_intersects_line(
            poly_ra, poly_dec, np.radians(line_ra), np.radians(line_dec)
        )
        == intersects
    )


def test_polygon_intersects_polygon():
    a = polygon_string_to_arrays("2:2, 2:3, 3:3, 3:2")
    b = polygon_string_to_arrays("1:1, 1:3, 3:3, 3:1")
    assert polygon_intersects_polygon(a[0], a[1], b[0], b[1])

    c = polygon_string_to_arrays("0:0, 0:1, 1:1, 1:0")
    assert not polygon_intersects_polygon(a[0], a[1], c[0], c[1])


def test_polygon_looped_around_polygon():
    a = polygon_string_to_arrays("2:2, 2:3, 3:3, 3:2")
    assert all(a[0] == np.radians((2, 2, 3, 3)))
    assert all(a[1] == np.radians((2, 3, 3, 2)))
    b = polygon_string_to_arrays("0:0, 0:5, 5:5, 5:0, 1:0, 1:1, 4:1, 4:4, 1:4, 1:0")
    assert not polygon_intersects_polygon(a[0], a[1], b[0], b[1])

    c = polygon_string_to_arrays("0:0, 0:5, 5:5, 5:0, 1:0, 1:2.5, 4:2.5, 4:4, 1:4, 1:0")
    assert polygon_intersects_polygon(a[0], a[1], c[0], c[1])


def test_polygon_intersects_cap():
    poly = polygon_string_to_arrays("2:2, 2:3, 3:3, 3:2")

    caps = [
        np.radians((0, 0, 1)),
        np.radians((0, 0, 3)),
        np.radians((2.5, 2.5, 0.2)),
        np.radians((2.5, 2.5, 1)),
    ]

    results = [
        (False, False, True, True),  # polygon contains center
        (False, False, True, False),  # polygon contains area
        (False, True, True, True),  # polygon intersects area
        (False, False, False, True),  # area contains polygon
    ]

    for intersection_type, intersection_results in zip(range(4), results):
        for i in range(4):
            result = polygon_intersects_cap(
                poly[0],
                poly[1],
                caps[i][0],
                caps[i][1],
                caps[i][2],
                intersection_type,
            )
            assert result == intersection_results[i]
