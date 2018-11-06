# Licensed with the 3-clause BSD license.  See LICENSE for details.
import sqlite3
import pytest
from logging import Logger

import numpy as np
from astropy.coordinates import SkyCoord

from .. import util
from ..db import SBDB
from ..config import Config

# tile 1/2 the sky
N_tiles = 10
ra_steps = np.linspace(0, 2 * np.pi, N_tiles + 1)
dec_steps = np.linspace(-np.pi / 2, np.pi / 2, N_tiles + 1)
sky_tiles = np.zeros((10, N_tiles**2))
for i in range(N_tiles):
    for j in range(N_tiles):
        sky_tiles[:, i * N_tiles + j] = (
            np.mean(ra_steps[i:i+2]),
            np.mean(dec_steps[j:j+2]),
            ra_steps[i], dec_steps[j],
            ra_steps[i], dec_steps[j+1],
            ra_steps[i+1], dec_steps[j+1],
            ra_steps[i+1], dec_steps[j])
del ra_steps, dec_steps, i, j


@pytest.fixture
def db():
    db = sqlite3.connect(':memory:', 5, 0, None, True, SBDB)
    db.verify_tables(Logger('test'))
    db.add_object('C/1995 O1')
    db.add_object('2P')

    obsids = range(N_tiles**2)
    start = 2458119.5 + np.arange(N_tiles**2) * 30 / 86400
    stop = start + 30 / 86400

    db.add_observations(columns=[obsids, start, stop] + list(sky_tiles))

    yield db
    db.close()


@pytest.fixture
def config():
    return Config()


class Test_SBDB:
    def test_verify_tables(self, db):
        logger = Logger('test')
        db.verify_tables(logger)
        c = db.execute('''
        SELECT count() FROM sqlite_master
        WHERE type='table'
        AND (
            name='obj' OR name='eph' OR name='eph_tree' OR
            name='obs' OR name='obs_tree' OR name='obs_found'
        )''')
        count = c.fetchone()[0]
        assert count == 6

        db.execute('drop table eph')
        db.verify_tables(logger)
        c = db.execute('''
        SELECT count() FROM sqlite_master
        WHERE type='table'
        AND (
            name='obj' OR name='eph' OR name='eph_tree' OR
            name='obs' OR name='obs_tree' OR name='obs_found'
        )''')
        count = c.fetchone()[0]
        assert count == 6

    def test_add_ephemeris_mpc_fixed(self, db):
        db.add_ephemeris(2, '500', 2458119.5, 2458121.5, step='1d',
                         source='mpc', cache=True)
        c = db.execute('select count() from eph').fetchone()[0]
        assert c == 3

    def test_add_ephemeris_mpc_variable(self, db):
        db.add_ephemeris(2, '500', 2457799.5, 2457809.5, step=None,
                         source='mpc', cache=True)
        c = db.execute('select count() from eph').fetchone()[0]
        assert c == 36

    def test_add_object(self, db):
        row = db.execute('select * from obj where desg="C/1995 O1"'
                         ).fetchone()
        assert row[0] == 1
        assert row[1] == 'C/1995 O1'

    def test_add_observations(self, db):
        c = db.execute('select count() from obs').fetchone()[0]
        assert c == N_tiles**2
        c = db.execute('select count() from obs_tree').fetchone()[0]
        assert c == N_tiles**2

    def test_get_ephemeris(self, db):
        db.add_ephemeris(2, '500', 2458119.5, 2458121.5, step='1d',
                         source='mpc', cache=True)
        eph = db.get_ephemeris(2, 2458119.5, 2458121.5)
        assert len(eph) == 3

    def test_get_ephemeris_exact(self, db):
        epochs = (2458119.5, 2458120.5, 2458121.5)
        eph = db.get_ephemeris_exact('2P', '500', epochs, source='jpl',
                                     cache=True)
        assert len(eph) == 3

    def test_get_ephemeris_interp(self, db):
        db.add_ephemeris(2, '500', 2458119.5, 2458121.5, step='1d',
                         source='mpc', cache=True)
        jdc = 2458120.0
        jda, jdb = 2458119.5, 2458120.5
        eph = db.get_ephemeris(2, jda, jdb)
        a = SkyCoord(eph[0]['ra'], eph[0]['dec'], unit='deg')
        b = SkyCoord(eph[1]['ra'], eph[1]['dec'], unit='deg')
        test = util.spherical_interpolation(a, b, jda, jdb, jdc)

        eph = db.get_ephemeris_interp(2, [jdc])
        assert np.isclose(eph.separation(test).value, 0)

    def test_get_ephemeris_segments(self, db):
        jda, jdb = 2458119.5, 2458121.5
        db.add_ephemeris(1, '500', jda, jdb, step='1d', source='mpc',
                         cache=True)
        db.add_ephemeris(2, '500', jda, jdb, step='1d', source='mpc',
                         cache=True)
        segments = db.get_ephemeris_segments()
        assert len(list(segments)) == 6

        segments = db.get_ephemeris_segments(start=jda, stop=jdb - 1)
        assert len(list(segments)) == 4

        segments = db.get_ephemeris_segments(objid=1, start=jda, stop=jdb - 1)
        assert len(list(segments)) == 2

    def test_clean_ephemeris(self, db):
        jda, jdb = 2458119.5, 2458121.5
        db.add_ephemeris(2, '500', jda, jdb, step='1d', source='mpc',
                         cache=True)
        eph = db.get_ephemeris(2, jda, jdb)
        assert len(eph) == 3
        count = db.clean_ephemeris(2, jda, jdb)
        assert count == 3
        eph = db.get_ephemeris(2, jda, jdb)
        assert len(eph) == 0

    def test_get_observations_by_id_errors(self, db):
        with pytest.raises(ValueError):
            db.get_observations_by_id(obsids=[1], columns='obsid')

    def test_get_observations_by_date(self, db):
        obsids = db.get_observations_by_date(2458119.5, 2458121.5,
                                             columns=['obsid'])
        assert len(obsids) == N_tiles**2

    def test_get_observations_by_id(self, db):
        obsids = db.get_observations_by_id(obsids=[1, 2, 3], generator=True)
        assert len(list(obsids)) == 3

    def test_get_observations_overlapping(self, db):
        db.add_ephemeris(2, '500', 2458119.5, 2458121.5, step='1d',
                         source='mpc', cache=True)
        eph = db.get_ephemeris(2, None, None)

        ra = [eph[i]['ra'] for i in range(len(eph))]
        dec = [eph[i]['dec'] for i in range(len(eph))]

        epochs = [eph[i]['jd'] for i in range(len(eph))]
        obs = db.get_observations_overlapping(ra=ra, dec=dec, epochs=epochs)

        # for N_tiles == 10, ephemeris will be in just one box
        assert len(obs) == 1

    def test_resolve_object(self, db):
        objid, desg = db.resolve_object(1)
        assert objid == 1
        assert desg == 'C/1995 O1'

        objid, desg = db.resolve_object('2P')
        assert objid == 2
        assert desg == '2P'
