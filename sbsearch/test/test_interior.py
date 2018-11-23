# Licensed with the 3-clause BSD license.  See LICENSE for details.
import pytest
import astropy.units as u
from ..util import RADec
from .. import interior


@pytest.mark.parametrize('point,test', (
    (RADec(0.5 * u.hourangle, 0.5 * u.deg), True),
    (RADec(-0.5 * u.hourangle, -0.5 * u.deg), False),
    (RADec(-0.5 * u.hourangle, 1.5 * u.deg), False),
    (RADec(0.5 * u.hourangle, 1.5 * u.deg), False),
    (RADec(0.5 * u.hourangle, -0.5 * u.deg), False)))
def test_interior_test(point, test):
    corners = RADec([0, 1, 1, 0] * u.hourangle, [0, 0, 1, 1] * u.deg)
    assert test == interior.interior_test(point, corners)

    corners = RADec([1, 1, 0, 0] * u.hourangle, [0, 1, 1, 0] * u.deg)
    assert test == interior.interior_test(point, corners)

    corners = RADec([0, 1, 0, 1] * u.hourangle, [0, 0, 1, 1] * u.deg)
    assert test == interior.interior_test(point, corners)
