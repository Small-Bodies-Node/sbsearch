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


@pytest.fixture
def sbs():
    with SBSearch(test=True) as sbs:
        sbs.clean_found(None)  # remove dummy found rows
        yield sbs


class TestSBSearch:
    def test_init(self):
        # just exercise
        with SBSearch() as sbs:
            assert sbs.db.engine is not None

    def test_disable_log(self, caplog):
        with SBSearch(disable_log=True, test=True) as sbs:
            level = sbs.logger.getEffectiveLevel()
            assert level == logging.ERROR

    def test_add_found(self, sbs):
        objid = sbs.db.resolve_object('2P')[0]
        obsids, foundids, tab = sbs.find_object(
            objid, vmax=25, save=False, update=False)
        obs = sbs.db.get_observations_by_id(obsids).all()
        foundids = sbs.add_found(objid, obs)
        assert len(obsids) == 1
        assert len(foundids) == 1
        assert obsids[0] in foundids

    def test_clean_ephemeris(self, sbs):
        eph = sbs.db.get_ephemeris(2, 2458119.5, 2458121.5).all()
        assert len(eph) == 3
        sbs.clean_ephemeris([2, '10P'])
        eph = sbs.db.get_ephemeris(2, 2458119.5, 2458121.5).all()
        assert len(eph) == 0

    def test_clean_found(self, sbs):
        foundids = sbs.add_found_by_id(2, [84], cache=True)

        count = sbs.clean_found([2, '10P'], 2458121.5, 2458122.5)
        assert count == 0

        found = sbs.db.get_found(2, 2458119.5, 2458122.5).all()
        assert len(found) == 1

        count = sbs.clean_found([2], 2458119.5, 2458122.5)
        assert count == 1

        found = sbs.db.get_found(2, 2458119.5, 2458122.5).all()
        assert len(found) == 0

    def test_clean_found_all(self, sbs):
        foundids = sbs.add_found_by_id(2, [84], cache=True)

        count = sbs.clean_found(None)
        assert count == 1

        found = sbs.db.get_found().all()
        assert len(found) == 0

    def test_find_by_ephemeris(self, sbs):
        eph = Ephem.from_dict({
            'ra': np.arange(10) * u.deg,
            'dec': np.arange(10) * u.deg,
            'Date': (2458119.5 + np.arange(10)) * u.day
        })
        found, tab = sbs.find_by_ephemeris(eph)
        assert len(found) == 4

    def test_find_by_ephemeris_single(self, sbs):
        eph = Ephem.from_dict({
            'ra': [1] * u.deg,
            'dec': [1] * u.deg,
            'Date': [2458119.5] * u.day
        })
        found, tab = sbs.find_by_ephemeris(eph)
        assert len(found) == 0

    def test_find_object(self, sbs):
        obsids, foundids, summary = sbs.find_object('2P', vmax=5)
        assert len(obsids) == 0  # too faint

        obsids, foundids, summary = sbs.find_object('2P', vmax=25)
        assert len(obsids) == 1

        obsids, foundids, summary = sbs.find_object(
            '2P', start=2458119.5, stop=2458121.5, vmax=25,
            source='jpl', cache=True)
        assert len(obsids) == 1

    def test_find_objects(self, sbs, caplog):
        # H-B has no ephemeris coverage, 2P should be found
        obsids, foundids, tab = sbs.find_objects(
            ['C/1995 O1', '2P'], save=True, cache=True)
        assert len(obsids) == 1

        # observation coverage, ephemeris coverage, but no found obseravtions
        obsids, foundids, tab = sbs.find_objects(
            ['2P'], start=2458119.53, stop=2458121.5, cache=True)
        assert tab is None

    def test_find_objects_messages(self, sbs, caplog):
        caplog.set_level(logging.INFO)
        sbs.logger = logging.getLogger('test dummy')
        sbs.find_objects(['2P'], cache=True)
        sbs.find_objects(['2P'], stop='2018-01-01', cache=True)
        sbs.find_objects(['2P'], start='2018-01-01', cache=True)
        sbs.find_objects(['2P'], start='2018-01-01', stop='2018-01-01',
                         cache=True)
        sbs.find_objects(['2P'], start='2018-01-01', stop='2018-01-02',
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

    def test_list_objects(self, sbs):
        tab = sbs.list_objects()
        assert all(tab['object ID'].data == [1, 2])
        assert all(tab['designation'].data == ['C/1995 O1', '2P'])

    def test_summarize_found(self, sbs):
        sbs.find_objects(['2P'], save=True, cache=True)
        tab = sbs.summarize_found()
        assert len(tab) == 1
        test = {
            'foundid': 4,
            'objid': 2,
            'obsid': 4,
            'jd': 2458119.529340278,
            'ra': 318.434,
            'dec': -17.707,
        }
        for k, v in test.items():
            if isinstance(v, int):
                assert tab[0][k] == v
            else:
                # coarse comparison in case ephemeris changes
                assert np.isclose(tab[0][k], v, rtol=0.01)

        tab = sbs.summarize_found(objects=['C/1995 O1'])
        assert tab is None

    def test_summarize_object_coverage(self, sbs):
        obsids, foundids, tab = sbs.find_objects(
            ['C/1995 O1', '2P'], save=True, cache=True)

        cov = sbs.summarize_object_coverage('eph', objects=[1, 2, '10P'],
                                            length=59)
        assert cov[0]['coverage'] == '-' * 59
        assert cov[2]['coverage'] == '-' * 59
        test = '+----------------------------+----------------------------+'
        assert cov[1]['coverage'] == test

        cov = sbs.summarize_object_coverage('found', length=60)
        assert cov[0]['coverage'] == '-' * 60
        test = '--+---------------------------------------------------------'
        assert cov[1]['coverage'] == test

    def test_summarize_object_coverage_error(self, sbs):
        with pytest.raises(ValueError):
            sbs.summarize_object_coverage('blah')

    def test_summarize_observations(self, sbs):
        tab = sbs.summarize_observations([987654321])
        assert tab is None

    def test_update_ephemeris(self, sbs):
        objid = sbs.db.resolve_object('2P')[0]
        start, stop = 2458119.5, 2458121.5
        N_eph = sbs.db.get_ephemeris(objid, None, None).count()
        assert N_eph == 3

        sbs.update_ephemeris(['2P'], start, stop, cache=True)
        N_eph = sbs.db.get_ephemeris(objid, None, None).count()
        assert N_eph == 3

        sbs.update_ephemeris(['10P'], start, stop, cache=True)
        objid = sbs.db.resolve_object('10P')[0]
        N_eph = sbs.db.get_ephemeris(objid, None, None).count()
        assert N_eph == 3
