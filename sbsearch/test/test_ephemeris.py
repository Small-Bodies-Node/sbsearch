# Licensed with the 3-clause BSD license.  See LICENSE for details.

from typing import List
import pytest

import numpy as np
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
from astropy.tests.helper import remote_data

from ..ephemeris import (get_ephemeris_generator, set_ephemeris_generator,
                         EphemerisGenerator, Horizons)
from ..target import MovingTarget
from ..model import Ephemeris


@remote_data
def test_get_set_ephemeris_generator():
    set_ephemeris_generator('jpl')
    g: EphemerisGenerator = get_ephemeris_generator()
    assert g == Horizons


@remote_data
def test_set_ephemeris_generator_error():
    with pytest.raises(ValueError):
        set_ephemeris_generator('invalid')


@remote_data
@pytest.fixture()
def encke() -> List[Ephemeris]:
    target: MovingTarget = MovingTarget('2P')

    g: Horizons = Horizons
    start: Time = Time('2020-06-01')
    stop: Time = Time('2020-07-01')
    step: u.Quantity = 10 * u.day
    return g.target_over_date_range(
        '500@', target, start, stop, step=step, cache=True)


@remote_data
class TestHorizons:
    def test_lt(self, encke):
        assert not (encke[0] < encke[0])
        assert encke[0] < encke[1]
        assert not (encke[1] < encke[0])

    def test_le(self, encke):
        assert encke[0] <= encke[0]
        assert encke[0] <= encke[1]
        assert not (encke[1] <= encke[0])

    def test_eq(self, encke):
        assert encke[0] == encke[0]
        assert not (encke[0] == encke[1])
        assert not (encke[1] == encke[0])

    def test_ne(self, encke):
        assert not (encke[0] != encke[0])
        assert encke[0] != encke[1]
        assert encke[1] != encke[0]

    def test_gt(self, encke):
        assert not (encke[0] > encke[0])
        assert not (encke[0] > encke[1])
        assert encke[1] > encke[0]

    def test_ge(self, encke):
        assert encke[0] >= encke[0]
        assert not (encke[0] >= encke[1])
        assert encke[1] >= encke[0]

    def test_at_dates(self):
        target: MovingTarget = MovingTarget('2P')
        # get these dates in reverse time order to be sure they are returned
        # in reverse order (Horizons requires ephemerides in order).
        dates: Time = Time(('2020-07-01', '2020-06-01'))
        g: Horizons = Horizons
        eph: List[Ephemeris] = g.target_at_dates(
            '500@', target, dates, cache=True)
        assert np.allclose([e.ra for e in eph], [120.97483, 65.4021],
                           rtol=1e-3)
        assert np.allclose([e.dec for e in eph], [17.72957, 26.36761],
                           rtol=1e-3)
        assert np.allclose([e.mjd for e in eph], dates.mjd)

    def test_at_dates_single(self):
        target: MovingTarget = MovingTarget('2P')
        date: Time = Time('1980-06-01')
        g: Horizons = Horizons
        eph: List[Ephemeris] = g.target_at_dates(
            '500@', target, date, cache=True)
        assert np.allclose([e.ra for e in eph], [18.94464], rtol=1e-3)
        assert np.allclose([e.dec for e in eph], [13.55878], rtol=1e-3)
        assert np.allclose([e.mjd for e in eph], date.mjd)

    def test_cometary_asteroid(self):
        target: MovingTarget = MovingTarget('174P')
        date: Time = Time('1980-06-01')
        g: Horizons = Horizons
        eph: List[Ephemeris] = g.target_at_dates(
            '500@', target, date, cache=True)
        assert len(eph) > 0

    def test_over_date_range(self, encke):
        assert np.allclose([e.ra for e in encke],
                           [65.4021, 81.36382, 100.97247, 120.97483],
                           rtol=1e-3)
        assert np.allclose([e.dec for e in encke],
                           [26.36761, 26.87021, 24.42504, 17.72957],
                           rtol=1e-3)
        assert np.allclose([e.mjd for e in encke],
                           [59001.0, 59011.0, 59021.0, 59031.0])

    @pytest.mark.parametrize(
        "start,stop",
        (
            ('2018-08-13', '2018-08-17'),
            ('2018-11-03', '2018-11-07')
        )
    )
    def test_over_date_range_adaptable(self, start, stop):
        """Wirtanen in 2018, closest approach 0.08 au.

        This test will work only with approaching objects.

        """
        target: MovingTarget = MovingTarget('46P')

        g: Horizons = Horizons
        start: Time = Time(start)
        stop: Time = Time(stop)
        eph: List[Ephemeris] = g.target_over_date_range('500@', target, start, stop,
                                                        cache=True)

        d: np.ndarray = np.diff([e.mjd for e in eph])
        delta: np.ndarray = np.array([e.delta for e in eph][1:])
        limit: float
        step: float
        for limit, step in zip([1, 0.25, 0], [1, 4 / 24, 1 / 24]):
            i: np.ndarray = delta > limit
            if i.sum() > 1:
                assert np.isclose(d[i][1:].mean(), step)
            if i.sum() > 0:
                delta = delta[~i]
                d = d[~i]
