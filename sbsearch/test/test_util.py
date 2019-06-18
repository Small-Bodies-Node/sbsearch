# Licensed with the 3-clause BSD license.  See LICENSE for details.
import struct
import sqlite3

import pytest
import numpy as np
import astropy.units as u
from sbpy.data import Ephem

from .. import util
from ..util import RADec, FieldOfView, Line, Point
from ..schema import Eph


class TestRADec:
    def test_init(self):
        ra = [0.477241085336078, 0.476386921200151, 0.475610424215939,
              0.47491560864072, 0.474305965132999]
        dec = [0.529339687106884, 0.525405714972889, 0.521444864768413,
               0.517461325283661, 0.513459808907613]
        a = RADec(ra, dec, unit='rad')
        b = RADec(a)
        assert np.allclose(b.ra.rad, ra)
        assert np.allclose(b.dec.rad, dec)

        c = RADec((ra[0], dec[0]), unit='rad')
        assert c.ra.rad == ra[0]
        assert c.dec.rad == dec[0]

    def test_init_error(self):
        with pytest.raises(ValueError):
            RADec('asdf', 'fdsa', 'somethingelse')

    def test_from_eph_single(self):
        eph = Eph(ra=1, dec=1)
        a = RADec.from_eph(eph)
        assert np.allclose((a.ra.deg, a.dec.deg), (1, 1))

    def test_from_eph_array(self):
        eph = (Eph(ra=1, dec=1), Eph(ra=2, dec=2))
        a = RADec.from_eph(eph)
        assert np.allclose(a.ra.deg, (1, 2))
        assert np.allclose(a.dec.deg, (1, 2))

    def test_repr(self):
        a = RADec(1, 0.5, unit='rad')
        assert repr(a) == '<RADec: ra=1.0 rad, dec=0.5 rad>'

    def test_len(self):
        a = RADec(1, 0.5, unit='rad')
        assert len(a) == 1

    def test_getitem(self):
        a = RADec([1], [1], unit='rad')
        b = a[0]
        assert np.allclose((a.ra[0].rad, a.dec[0].rad), (b.ra.rad, b.dec.rad))

    def test_separation(self):
        a = RADec([1], [0.5], unit='rad')
        b = RADec([1.01], [0.52], unit='rad')
        assert np.isclose(a.separation(b).to('rad').value,
                          0.021821164590249686)

    def test_xyz(self):
        from numpy import pi
        a = RADec([0, pi/2, pi, 0], [0, 0, 0, pi/2], unit='rad')
        assert np.allclose(a.xyz, ((1, 0, -1, 0), (0, 1, 0, 0), (0, 0, 0, 1)))


class TestFieldOfView:
    def test_init_error(self):
        with pytest.raises(TypeError):
            FieldOfView([1, 2])

    def test_str(self):
        coords = RADec(((0, 0), (0, 1), (1, 1), (1, 0)), unit='deg')
        fov = str(FieldOfView(coords))
        assert fov == 'SRID=40001;POLYGON((0.0 0.0,0.0 1.0,1.0 1.0,1.0 0.0,0.0 0.0))'


class TestLine:
    def test_init_error(self):
        with pytest.raises(TypeError):
            Line([1, 2])

    def test_str(self):
        coords = RADec(((0, 0), (0, 1), (1, 1), (1, 0)), unit='deg')
        line = str(Line(coords))
        assert line == 'SRID=40001;LINESTRING((0.0 0.0,0.0 1.0,1.0 1.0,1.0 0.0))'

    def test_from_eph(self):
        eph = [Eph(ra=1, dec=1), Eph(ra=2, dec=2)]
        line = str(Line.from_eph(eph))
        assert line == 'SRID=40001;LINESTRING((1.0 1.0,2.0 2.0))'

    def test_from_ephem(self):
        eph = Ephem(dict(ra=[1, 2] * u.deg, dec=[1, 2] * u.deg))
        line = str(Line.from_ephem(eph))
        assert line == 'SRID=40001;LINESTRING((1.0 1.0,2.0 2.0))'


class TestPoint:
    def test_init_error(self):
        with pytest.raises(TypeError):
            Point([1, 2])

    def test_str(self):
        coords = RADec((0, 0), unit='deg')
        point = str(Point(coords))
        assert point == 'SRID=40001;POINT(0.0 0.0)'

    def test_from_eph(self):
        eph = Eph(ra=1, dec=1)
        point = str(Point.from_eph(eph))
        assert point == 'SRID=40001;POINT(1.0 1.0)'

    def test_from_ephem(self):
        eph = Ephem(dict(ra=[1] * u.deg, dec=[1] * u.deg))
        point = str(Point.from_ephem(eph))
        assert point == 'SRID=40001;POINT(1.0 1.0)'


def test_epochs_to_time():
    t = util.epochs_to_time(['2018-01-01', 2455000.5])
    assert np.allclose(t.jd, (2458119.5, 2455000.5))


def test_epochs_to_jd():
    jd = util.epochs_to_jd(['2018-01-01', 2455000.5])
    assert np.allclose(jd, (2458119.5, 2455000.5))


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
    assert np.isclose(c2.ra.rad, 0)
    assert np.isclose(c2.dec.rad, 0)

    c2 = util.spherical_interpolation(c0, c1, 0, 2, 0)
    assert np.isclose(c2.ra.rad, c0.ra.rad)
    assert np.isclose(c2.dec.rad, c0.dec.rad)

    c2 = util.spherical_interpolation(c0, c1, 0, 2, 2)
    assert np.isclose(c2.ra.rad, c1.ra.rad)
    assert np.isclose(c2.dec.rad, c1.dec.rad)

    c2 = util.spherical_interpolation(c0, c1, 0, 2, 0)
    assert np.isclose(c2.ra.rad, c0.ra.rad)
    assert np.isclose(c2.dec.rad, c0.dec.rad)

    c2 = util.spherical_interpolation(c0, c1, 0, 2, 2)
    assert np.isclose(c2.ra.rad, c1.ra.rad)
    assert np.isclose(c2.dec.rad, c1.dec.rad)


def test_spherical_interpolation_error():
    c0 = RADec(-0.1, 0.1, unit='rad')
    c1 = RADec(0.1, -0.1, unit='rad')
    with pytest.raises(ValueError):
        util.spherical_interpolation(c0, c1, 0, 0, 0)


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
