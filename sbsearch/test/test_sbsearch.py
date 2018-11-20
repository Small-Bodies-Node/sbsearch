# Licensed with the 3-clause BSD license.  See LICENSE for details.
import sqlite3
import pytest
from logging import Logger
from itertools import repeat

import numpy as np
import astropy.time as Time

from .. import util
from ..sbsearch import SBSearch
from ..config import Config
from .skytiles import sky_tiles, N_tiles


@pytest.fixture
def sbs():
    config = Config(database=':memory:')
    with SBSearch(config) as sbs:
        sbs.db.add_object('C/1995 O1')
        objid = sbs.db.add_object('2P')

        obsids = range(N_tiles**2)
        start = 2458119.5 + np.arange(N_tiles**2) * 30 / 86400
        stop = start + 30 / 86400

        columns = [obsids, repeat('test'), start, stop, sky_tiles]
        sbs.db.add_observations(zip(*columns))
        sbs.update_ephemeris([objid], 2458119.5, 2458121.5, step='1d',
                             cache=True)

        yield sbs


class TestSBSearch:
    def test_update_ephemeris(self, sbs):
        objid = sbs.db.resolve_object('2P')[0]
        start, stop = 2458119.5, 2458121.5
        N_eph = len(sbs.db.get_ephemeris(objid, None, None))
        N_eph_tree = len(list(sbs.db.get_ephemeris_segments(
            objid=objid, start=None, stop=None)))
        assert N_eph == 3
        assert N_eph_tree == 3

        sbs.update_ephemeris([objid], start, stop, cache=True)
        N_eph = len(sbs.db.get_ephemeris(objid, None, None))
        N_eph_tree = len(list(sbs.db.get_ephemeris_segments(
            objid=objid, start=None, stop=None)))
        assert N_eph == 3
        assert N_eph_tree == 3

    def test_find_object(self, sbs):
        obsids = sbs.find_object('2P', vmax=5)
        assert len(obsids) == 0  # too faint

        obsids = sbs.find_object('2P', vmax=25)
        assert len(obsids) == 1
