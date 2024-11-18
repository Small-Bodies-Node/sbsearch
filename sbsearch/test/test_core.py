# Licensed with the 3-clause BSD license.  See LICENSE for details.

import pytest

from typing import List
import numpy as np

from .. import core
from ..model.core import Ephemeris, Observation
from ..exceptions import EdgeCrossingDetected


def test_polygon_string_to_arrays():
    ra, dec = core.polygon_string_to_arrays("1:0, 10:20, 101:0")
    assert np.all(ra == np.radians((1, 10, 101)))
    assert np.all(dec == np.radians((0, 20, 0)))


def test_test_edges():
    # no error raised
    ra = [0, 0, 1, 1]
    dec = [0, 1, 1, 0]
    obs = Observation()

    # set FOV also runs verify_corner_order
    obs.set_fov(ra, dec)

    # reverse vertices, should still be good
    obs.set_fov(ra[::-1], dec[::-1])

    # change the order and now we should get an exception
    ra = [0, 0, 1, 1]
    dec = [0, 1, 0, 1]
    with pytest.raises(EdgeCrossingDetected):
        obs.set_fov(ra, dec)

    # this one is actually a line, it might work with S2, but we don't want it
    # in SBSearch
    ra = [0, 1, 1, 0]
    dec = [0, 1, 1, 0]
    with pytest.raises(EdgeCrossingDetected):
        obs.set_fov(ra, dec)

    # FOVs covering a pole should not be a problem

    # this one is OK
    ra = [0, 90, 180, 270]
    dec = [88, 88, 88, 88]
    obs.set_fov(ra, dec)

    # this one is not
    ra = [0, 90, 270, 180]
    with pytest.raises(EdgeCrossingDetected):
        obs.set_fov(ra, dec)


def test_ephemeris_uncertainty_offsets():
    eph: List[Ephemeris] = [
        Ephemeris(ra=0, dec=0, dra=0, ddec=1e-3, unc_a=30, unc_b=10, unc_theta=0)
    ]

    with pytest.raises(ValueError):
        core.ephemeris_uncertainty_offsets(eph)

    eph.append(
        Ephemeris(ra=0, dec=1, dra=0, ddec=1e-3, unc_a=30, unc_b=10, unc_theta=0)
    )
    a: np.ndarray  # parallel to target motion
    b: np.ndarray  # perpendicular to target motion
    a, b = core.ephemeris_uncertainty_offsets(eph)

    # unc_a is 30 along PA=0, unc_a --> parallel --> a
    assert np.allclose(a, np.radians(30 / 3600))
    assert np.allclose(b, np.radians(10 / 3600))

    # rotate unc ellipse
    for e in eph:
        e.unc_theta = 45

    # now unc_a dominates:
    a, b = core.ephemeris_uncertainty_offsets(eph)
    assert np.allclose(a, np.radians(30 * np.cos(np.radians(45)) / 3600))
    assert np.allclose(b, np.radians(30 * np.cos(np.radians(45)) / 3600))

    # rotate unc ellipse
    for e in eph:
        e.unc_theta = 90

    # unc_a is 30 along PA=90, unc_a --> perpendicular --> b
    a, b = core.ephemeris_uncertainty_offsets(eph)
    assert np.allclose(a, np.radians(10 / 3600))
    assert np.allclose(b, np.radians(30 / 3600))
