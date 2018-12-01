# Licensed with the 3-clause BSD license.  See LICENSE for details.
import sqlite3
import pytest
import logging
from itertools import repeat

import numpy as np
import astropy.units as u
import astropy.time as Time
from sbpy.data import Ephem

from .. import util
from ..sbsearch import SBSearch
from ..config import Config
from .skytiles import sky_tiles, N_tiles


def sbs_for_testing():
    config = Config(database=':memory:')
    sbs = SBSearch(config)
    sbs.db.add_object('C/1995 O1')
    objid = sbs.db.add_object('2P')

    obsids = range(N_tiles**2)
    start = 2458119.5 + np.arange(N_tiles**2) * 30 / 86400
    stop = start + 30 / 86400

    columns = [obsids, repeat('test'), start, stop, sky_tiles]
    sbs.db.add_observations(zip(*columns))
    sbs.update_ephemeris([objid], 2458119.5, 2458121.5, step='1d',
                         cache=True)

    return sbs


@pytest.fixture
def sbs():
    sbs = sbs_for_testing()
    yield sbs
    sbs.__exit__()


class TestSBSearch:
    def test_disable_log(self, caplog):
        config = Config(database=':memory:')
        with SBSearch(config, disable_log=True) as sbs:
            level = sbs.logger.getEffectiveLevel()
            assert level == logging.ERROR

    def test_add_found(self, sbs):
        objid = sbs.db.resolve_object('2P')[0]
        obs = sbs.find_object(objid, vmax=25)
        foundids = sbs.add_found(objid, obs, '500', cache=True)
        assert len(foundids) == 1
        assert foundids[0] == 1

    def test_available_objects(self, sbs):
        tab = sbs.available_objects()
        assert all(tab['object ID'].data == [1, 2])
        assert all(tab['designation'].data == ['C/1995 O1', '2P'])

    def test_clean_ephemeris(self, sbs):
        eph = sbs.db.get_ephemeris(2, 2458119.5, 2458121.5)
        assert len(eph) == 3
        sbs.clean_ephemeris([2, '10P'])
        eph = sbs.db.get_ephemeris(2, 2458119.5, 2458121.5)
        assert len(eph) == 0

    def test_clean_found(self, sbs):
        foundids = sbs.add_found_by_id(2, [84], '500', cache=True)
        sbs.clean_found([2, '10P'], 2458121.5, 2458122.5)
        found = sbs.db.get_found(2, 2458119.5, 2458122.5)
        assert len(found) == 1
        sbs.clean_found([2], 2458119.5, 2458122.5)
        found = sbs.db.get_found(2, 2458119.5, 2458122.5)
        assert len(found) == 0

    def test_clean_found_all(self, sbs):
        foundids = sbs.add_found_by_id(2, [84], '500', cache=True)
        sbs.clean_found(None)
        found = sbs.db.get_found()
        assert len(found) == 0

    def test_find_by_ephemeris(self, sbs):
        eph = Ephem.from_dict({
            'ra': np.arange(10) * u.deg,
            'dec': np.arange(10) * u.deg,
            'Date': (2458119.5 + np.arange(10)) * u.day
        })
        n, found, tab = sbs.find_by_ephemeris(eph)
        assert len(found) == 1

    def test_find_by_ephemeris_error(self, sbs):
        eph = Ephem.from_dict({
            'ra': np.arange(2) * u.deg,
            'dec': np.arange(2) * u.deg,
            'Date': (2458119.5 + np.arange(2)) * u.day
        })
        with pytest.raises(ValueError):
            n, found, tab = sbs.find_by_ephemeris(eph)

    def test_find_object(self, sbs):
        obs = sbs.find_object('2P', vmax=5)
        assert len(obs) == 0  # too faint

        obs = sbs.find_object('2P', vmax=25)
        assert len(obs) == 1

    def test_find_objects(self, sbs, caplog):
        # H-B has no ephemeris coverage, 2P should be found
        tab = sbs.find_objects(['C/1995 O1', '2P'], save=True, cache=True)
        assert len(tab) == 1

        # observation coverage, ephemeris coverage, but no found obseravtions
        tab = sbs.find_objects(['2P'], start=2458119.53, stop=2458121.5,
                               cache=True)
        assert len(tab) == 0

    def test_find_objects_messages(self, sbs, caplog):
        caplog.set_level(logging.INFO)
        sbs.logger = logging.getLogger('test dummy')
        tab = sbs.find_objects(['2P'], cache=True)
        tab = sbs.find_objects(['2P'], stop='2018-01-01', cache=True)
        tab = sbs.find_objects(['2P'], start='2018-01-01', cache=True)
        tab = sbs.find_objects(['2P'], start='2018-01-01', stop='2018-01-01',
                               cache=True)
        tab = sbs.find_objects(['2P'], start='2018-01-01', stop='2018-01-02',
                               cache=True)
        expected = (
            'in all observations',
            'in observations ending 2018-01-01 00:00 UT',
            'in observations starting 2018-01-01 00:00 UT',
            'in observations on 2018-01-01',
            'in observations from 2018-01-01 00:00 UT to 2018-01-02 00:00 UT'
        )
        captured = '\n'.join([r.getMessage() for r in caplog.records])
        for test in expected:
            assert test in captured

    def test_found_summary(self, sbs):
        sbs.find_objects(['2P'], save=True, cache=True)
        tab = sbs.found_summary()
        assert len(tab) == 1
        row = tuple(tab[0])
        test = (1, '2P', 2458119.529340278, 318.434, -17.707,
                0.104, 0.086, 20.53, 3.36, 8.900, 4.11, 9.7, 35.2)
        assert row[0] == test[0]
        assert row[1] == test[1]
        # coarse comparison in case ephemeris changes
        assert np.allclose(row[2:], test[2:], rtol=0.01)

        tab = sbs.found_summary(objects=['C/1995 O1'])
        assert tab is None

    def test_object_coverage(self, sbs):
        tab = sbs.find_objects(['C/1995 O1', '2P'], save=True, cache=True)

        cov = sbs.object_coverage('eph', objects=[1, 2, '10P'], length=59)
        assert cov[0]['coverage'] == '-' * 59
        assert cov[2]['coverage'] == '-' * 59
        test = '+----------------------------+----------------------------+'
        assert cov[1]['coverage'] == test

        cov = sbs.object_coverage('found', length=60)
        print(cov.meta)
        assert cov[0]['coverage'] == '-' * 60
        test = '--------------------------------------------------+---------'
        assert cov[1]['coverage'] == test

    def test_object_coverage_error(self, sbs):
        with pytest.raises(ValueError):
            sbs.object_coverage('blah')

    def test_update_ephemeris(self, sbs):
        objid = sbs.db.resolve_object('2P')[0]
        start, stop = 2458119.5, 2458121.5
        N_eph = len(sbs.db.get_ephemeris(objid, None, None))
        ephids, segments = sbs.db.get_ephemeris_segments(
            objid=objid, start=None, stop=None)
        N_eph_tree = len(ephids)
        assert N_eph == 3
        assert N_eph_tree == 3

        sbs.update_ephemeris(['2P'], start, stop, cache=True)
        N_eph = len(sbs.db.get_ephemeris(objid, None, None))
        ephids, segments = sbs.db.get_ephemeris_segments(
            objid=objid, start=None, stop=None)
        N_eph_tree = len(ephids)
        assert N_eph == 3
        assert N_eph_tree == 3

        sbs.update_ephemeris(['10P'], start, stop, cache=True)
        objid = sbs.db.resolve_object('10P')[0]
        N_eph = len(sbs.db.get_ephemeris(objid, None, None))
        ephids, segments = sbs.db.get_ephemeris_segments(
            objid=objid, start=None, stop=None)
        N_eph_tree = len(ephids)
        assert N_eph == 3
        assert N_eph_tree == 3
