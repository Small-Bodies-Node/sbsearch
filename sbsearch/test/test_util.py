# Licensed with the 3-clause BSD license.  See LICENSE for details.
import sqlite3

import pytest
import numpy as np
import astropy.units as u
from sbpy.data import Ephem

from .. import util
from ..util import RADec


def test_assemble_sql():
    cmd = 'SELECT * FROM table'
    parameters = []
    constraints = [('v > ?', 1), ('v < 2', None)]
    r = util.assemble_sql(cmd, parameters, constraints)
    assert r[0] == 'SELECT * FROM table WHERE v > ? AND v < 2'
    assert r[1] == [1]


def test_eph_to_limits():
    from numpy import pi
    ra = [0, 0, 0]
    dec = [-pi / 2, 0, pi / 2]
    eph = RADec(ra, dec, unit='rad')
    jd = [2400000.5, 2400001.5, 2400002.5]
    half_step = 0.5 * u.day
    r = util.eph_to_limits(eph, jd, half_step)
    assert np.allclose(r, [0.5, 1.5, 0.707107, 1, 0, 0, -0.707107, 0.707107])


def test_eph_to_limits_continuity():
    ra = [0.477241085336078, 0.476386921200151, 0.475610424215939,
          0.47491560864072, 0.474305965132999]
    dec = [0.529339687106884, 0.525405714972889, 0.521444864768413,
           0.517461325283661, 0.513459808907613]
    jd = [2458429.5, 2458430.5, 2458431.5, 2458432.5, 2458433.5]
    coords = RADec(ra, dec, unit='rad')
    half_step = 0.5 * u.day
    limits = []
    for i in [1, 2, 3]:
        if i == 0:
            indices = (1, 0, 1)
        elif i == len(coords) - 1:
            indices = (-2, -1, -2)
        else:
            indices = (i - 1, i, i + 1)

        print(indices)
        c = tuple((coords[j] for j in indices))
        _jd = tuple((jd[j] for j in indices))
        limits.append(util.eph_to_limits(c, _jd, half_step))
        print(limits[-1])

    # mjd, x, and y edges should match (z will not)
    assert np.isclose(limits[0][1], limits[1][0])
    assert np.isclose(limits[1][1], limits[2][0])
    assert np.isclose(limits[0][3], limits[1][2])
    assert np.isclose(limits[1][3], limits[2][2])
    assert np.isclose(limits[0][5], limits[1][4])
    assert np.isclose(limits[1][5], limits[2][4])


def test_epochs_to_time():
    t = util.epochs_to_time(['2018-01-01', 2455000.5])
    assert np.allclose(t.jd, (2458119.5, 2455000.5))


def test_epochs_to_jd():
    jd = util.epochs_to_jd(['2018-01-01', 2455000.5])
    assert np.allclose(jd, (2458119.5, 2455000.5))


def test_date_constraints():
    constraints = util.date_constraints(1, 2)
    assert constraints == [('jd>=?', 1), ('jd<=?', 2)]


@pytest.mark.parametrize('point,test', (
    (RADec(0.5 * u.hourangle, 0.5 * u.deg), True),
    (RADec(-0.5 * u.hourangle, -0.5 * u.deg), False),
    (RADec(-0.5 * u.hourangle, 1.5 * u.deg), False),
    (RADec(0.5 * u.hourangle, 1.5 * u.deg), False),
    (RADec(0.5 * u.hourangle, -0.5 * u.deg), False)))
def test_interior_test(point, test):
    corners = RADec([0, 1, 1, 0] * u.hourangle, [0, 0, 1, 1] * u.deg)
    assert test == util.interior_test(point, corners)

    corners = RADec([1, 1, 0, 0] * u.hourangle, [0, 1, 1, 0] * u.deg)
    assert test == util.interior_test(point, corners)

    corners = RADec([0, 1, 0, 1] * u.hourangle, [0, 0, 1, 1] * u.deg)
    assert test == util.interior_test(point, corners)


def test_iterate_over():
    db = sqlite3.connect(':memory:')
    db.execute('CREATE TABLE t(a,b,c)')
    N = 10000
    db.executemany('INSERT INTO t VALUES (?,?,?)',
                   np.random.rand(N * 3).reshape(N, 3))
    c = db.execute('SELECT * FROM t')
    count = 0
    for row in util.iterate_over(c):
        count += 1

    assert count == N


def test_rd2xyz():
    from numpy import pi
    ra = [0, pi / 2, pi, 3 * pi / 2, 0, 0]
    dec = [0, 0, 0, 0, pi / 2, -pi / 2]
    xyz = util.rd2xyz(ra, dec)
    test = np.array(((1, 0, -1, 0, 0, 0),
                     (0, 1, 0, -1, 0, 0),
                     (0, 0, 0, 0, 1, -1)))
    assert np.allclose(xyz, test)


def test_spherical_interpolation():
    c0 = RADec(-0.1, 0.1, unit='rad')
    c1 = RADec(0.1, -0.1, unit='rad')
    c2 = util.spherical_interpolation(c0, c1, 0, 2, 1)
    assert np.isclose(c2.ra, 0)
    assert np.isclose(c2.dec, 0)

    c2 = util.spherical_interpolation(c0, c1, 0, 2, 0)
    assert np.isclose(c2.ra, c0.ra)
    assert np.isclose(c2.dec, c0.dec)

    c2 = util.spherical_interpolation(c0, c1, 0, 2, 2)
    assert np.isclose(c2.ra, c1.ra)
    assert np.isclose(c2.dec, c1.dec)


def test_vector_rotate():
    a = np.r_[1, 0, 0]
    n = np.r_[0, 0, 1]
    b = util.vector_rotate(a, n, np.pi / 2)
    assert np.allclose(b, [0, 1, 0])


def test_vmag_from_eph():
    eph = Ephem.from_dict({
        'Tmag': np.zeros(3) * u.mag,
        'Nmag': np.zeros(3) * u.mag,
        'V': np.zeros(3) * u.mag
    })

    v = util.vmag_from_eph(eph)
    assert np.allclose(v, 99)

    v = util.vmag_from_eph(eph, ignore_zero=False)
    assert np.allclose(v, 0)

    eph['Tmag'][0] = 1 * u.mag
    eph['Nmag'][0] = 2 * u.mag
    eph['V'][0] = 3 * u.mag

    eph['Nmag'][1] = 2 * u.mag
    eph['V'][1] = 3 * u.mag

    eph['V'][2] = 3 * u.mag
    v = util.vmag_from_eph(eph)
    assert np.allclose(v, [1, 2, 3])
