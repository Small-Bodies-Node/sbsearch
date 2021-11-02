# Licensed with the 3-clause BSD license.  See LICENSE for details.

import pytest

from typing import List
import numpy as np

from .. import core
from .. model.core import Ephemeris


def test_ephemeris_uncertainty_offsets():
    eph: List[Ephemeris] = [
        Ephemeris(ra=0, dec=0, dra=0, ddec=1e-3,
                  unc_a=30, unc_b=10, unc_theta=0)
    ]

    with pytest.raises(ValueError):
        core.ephemeris_uncertainty_offsets(eph)

    eph.append(Ephemeris(ra=0, dec=1, dra=0, ddec=1e-3,
                         unc_a=30, unc_b=10, unc_theta=0))
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
