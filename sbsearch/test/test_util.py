# Licensed with the 3-clause BSD license.  See LICENSE for details.

import sqlite3

import numpy as np
from astropy.coordinates import SkyCoord

from .. import util


def test_assemble_sql():
    cmd = 'SELECT * FROM table'
    parameters = []
    constraints = [('v > ?', 1), ('v < 2', None)]
    r = util.assemble_sql(cmd, parameters, constraints)
    assert r[0] == 'SELECT * FROM table WHERE v > ? AND v < 2'
    assert r[1] == [1]


def test_date_constraints():
    constraints = util.date_constraints(1, 2)
    assert constraints == [('jd>=?', 1), ('jd<=?', 2)]


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


def test_interior_test():
    pass


def test_rd2xyz():
    from numpy import pi
    ra = [0, pi / 2, pi, 3 * pi / 2, 0, 0]
    dec = [0, 0, 0, 0, pi / 2, -pi / 2]
    xyz = util.rd2xyz(ra, dec)
    test = np.array(((1, 0, -1, 0, 0, 0),
                     (0, 1, 0, -1, 0, 0),
                     (0, 0, 0, 0, 1, -1)))
    assert np.allclose(xyz, test)


def test_sc2xyz():
    from numpy import pi
    ra = [0, pi / 2, pi, 3 * pi / 2, 0, 0]
    dec = [0, 0, 0, 0, pi / 2, -pi / 2]
    c = SkyCoord(ra, dec, unit='rad')

    xyz = util.sc2xyz(c)
    test = np.array(((1, 0, -1, 0, 0, 0),
                     (0, 1, 0, -1, 0, 0),
                     (0, 0, 0, 0, 1, -1)))
    assert np.allclose(xyz, test)


def test_spherical_interpolation():
    c0 = SkyCoord(-0.1, 0.1, unit='rad')
    c1 = SkyCoord(0.1, -0.1, unit='rad')
    c2 = util.spherical_interpolation(c0, c1, 0, 2, 1)
    assert np.isclose(c2.ra.value, 0)
    assert np.isclose(c2.dec.value, 0)
