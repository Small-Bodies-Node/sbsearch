# Licensed with the 3-clause BSD license.  See LICENSE for details.
import pytest

import numpy as np
import astropy.units as u
from astropy.time import Time
from sbpy.data import Orbit

from .. import ephem


def test_format_epochs_error():
    with pytest.raises(ValueError):
        ephem._format_epochs((2, 1))


def test_generate():
    epochs = (2458119.5, 2458120.5, 2458121.5)
    eph = ephem.generate('2P', '500', epochs, source='jpl', cache=True)
    assert len(eph) == 3

    # this one has ephemeris uncertainties
    eph = ephem.generate('2019 AQ9', '500', epochs, source='mpc',
                         cache=True)
    assert len(eph) == 3

    # Encke does not
    eph = ephem.generate('2P', '500', epochs, source='mpc', cache=True)
    assert len(eph) == 3


def test_generate_range_fixed_steps():
    epochs = {
        'start': Time(2458200.5, format='jd').iso,
        'stop': Time(2458210.5, format='jd').iso,
        'step': '1d'
    }
    eph = ephem.generate('2P', '500', epochs, source='jpl', cache=True)
    assert len(eph) == 11


def test_generate_range_adaptable_steps():
    for source in ('jpl', 'mpc'):
        epochs = {
            'start': Time(2457799.5, format='jd').iso,
            'stop': Time(2457809.5, format='jd').iso,
            'step': None
        }
        eph = ephem.generate('2P', '500', epochs, source=source, cache=True)
        assert len(eph) == 36


def test_generate_from_orbit():
    epochs = (2458119.5, 2458120.5, 2458121.5)
    orbit = Orbit.from_dict({
        'targetname': ['2P'],
        'a': [2.215134573264697] * u.au,
        'e': [0.8483251746071773],
        'i': [11.78183005207527] * u.deg,
        'Omega': [334.5678392905074] * u.deg,
        'w': [186.5437009154704] * u.deg,
        'M': [143.2720471022976] * u.deg,
        'epoch': [2457097.5],
        'timescale': ['UTC'],
        'H': [14.2],
        'G': [0.15]
    })
    eph = ephem.generate('2P', '500', epochs, source='oorb', orbit=orbit)
    assert len(eph) == 3


def test_generate_error():
    epochs = (2458119.5, 2458120.5, 2458121.5)
    with pytest.raises(ValueError):
        ephem.generate('2P', '500', epochs, source='horizons')

    epochs = (2458119.5, 2458119.5, 2458119.5)
    with pytest.raises(ValueError):
        ephem.generate('2P', '500', epochs, source='jpl')


def test_generate_365():
    ephem.STEP_LIMIT = 300
    epochs = np.arange(365) + 2458119.5
    eph = ephem.generate('2P', '500', epochs, source='jpl', cache=True)
    assert len(eph) == 365


def test_generate_orbit():
    orb = ephem.generate_orbit('2P', [2458119.5], cache=True)
    assert len(orb) == 1

    epochs = {'start': '2018-01-01', 'stop': '2018-01-03', 'step': '1d'}
    orb = ephem.generate_orbit('2P', epochs, cache=True)
    assert len(orb) == 3


def test_generate_orbit_365():
    ephem.STEP_LIMIT = 300
    epochs = 2458119.5 + np.arange(365)
    orb = ephem.generate_orbit('2P', epochs, cache=True)
    assert len(orb) == 365


def test_generate_orbit_error():
    epochs = [2458119.5] * 3
    with pytest.raises(ValueError):
        orb = ephem.generate_orbit('2P', epochs, cache=True)
