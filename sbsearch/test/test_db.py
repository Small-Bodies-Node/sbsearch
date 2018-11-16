# Licensed with the 3-clause BSD license.  See LICENSE for details.
import sqlite3
import pytest
from logging import Logger
from itertools import repeat

import numpy as np
import astropy.time as Time

from .. import util
from ..util import RADec
from ..db import SBDB
from ..exceptions import BadObjectID

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
    db.verify_database(Logger('test'))
    db.add_object('C/1995 O1')
    db.add_object('2P')

    obsids = range(N_tiles**2)
    start = 2458119.5 + np.arange(N_tiles**2) * 30 / 86400
    stop = start + 30 / 86400

    columns = [obsids, repeat('test'), obsids, start, stop] + list(sky_tiles)
    db.add_observations(columns=columns)
    db.add_ephemeris(2, '500', 2458119.5, 2458121.5, step='1d',
                     source='mpc', cache=True)
    db.add_found(2, [1, 2, 3], '500', cache=True)

    yield db
    db.close()


class Test_SBDB:
    def test_verify_database(self, db):
        logger = Logger('test')
        db.verify_database(logger)
        c = db.execute('''
        SELECT count() FROM sqlite_master
        WHERE type='table'
        AND (
            name='obj' OR name='eph' OR name='eph_tree' OR
            name='obs' OR name='obs_tree' OR name='found'
        )''')
        count = c.fetchone()[0]
        assert count == 6

        db.execute('drop table eph')
        db.verify_database(logger)
        c = db.execute('''
        SELECT count() FROM sqlite_master
        WHERE type='table'
        AND (
            name='obj' OR name='eph' OR name='eph_tree' OR
            name='obs' OR name='obs_tree' OR name='found'
        )''')
        count = c.fetchone()[0]
        assert count == 6

    def test_make_test_db(self):
        # just exercise the code
        db = SBDB.make_test_db()
        assert db is not None

    def test_add_ephemeris_mpc_fixed(self, db):
        c = db.execute('select count() from eph').fetchone()[0]
        assert c == 3

    def test_add_ephemeris_mpc_variable(self, db):
        db.add_ephemeris(2, '500', 2457799.5, 2457809.5, step=None,
                         source='mpc', cache=True)
        c = db.execute('select count() from eph').fetchone()[0]
        assert c == 39  # 36 here + 3 add at top

    def test_add_found(self, db):
        rows = db.execute('select * from found').fetchall()
        assert len(rows) == 3

    def test_add_object(self, db):
        row = db.execute('select * from obj where desg="C/1995 O1"'
                         ).fetchone()
        assert row[0] == 1
        assert row[1] == 'C/1995 O1'

    def test_add_object_error(self, db):
        with pytest.raises(ValueError):
            db.add_object(1)

    def test_add_observations(self, db):
        c = db.execute('select count() from obs').fetchone()[0]
        assert c == N_tiles**2
        c = db.execute('select count() from obs_tree').fetchone()[0]
        assert c == N_tiles**2

        # add rows
        rows = [[100, 'test', 101, 5, 10, 0, 0, 1, 1, 1, -1, -1, -1, -1, 1]]
        db.add_observations(rows=rows)
        c = db.execute('select count() from obs').fetchone()[0]
        assert c == N_tiles**2 + 1

    def test_clean_ephemeris(self, db):
        jda, jdb = 2458119.5, 2458121.5
        eph = db.get_ephemeris(2, jda, jdb)
        assert len(eph) == 3
        count = db.clean_ephemeris(2, jda, jdb)
        assert count == 3
        eph = db.get_ephemeris(2, jda, jdb)
        assert len(eph) == 0

    def test_clean_found(self, db):
        count = db.clean_found(2, None, None)
        assert count == 3

    def test_get_ephemeris(self, db):
        eph = db.get_ephemeris(2, 2458119.5, 2458121.5)
        assert len(eph) == 3

    def test_get_ephemeris_all(self, db):
        eph = db.get_ephemeris(None, None, None, generator=True)
        assert len(list(eph)) == 3

    def test_get_ephemeris_exact(self, db):
        epochs = (2458119.5, 2458120.5, 2458121.5)
        eph = db.get_ephemeris_exact('2P', '500', epochs, source='jpl',
                                     cache=True)
        assert len(eph) == 3

    def test_get_ephemeris_exact_error(self, db):
        epochs = (2458119.5, 2458120.5, 2458121.5)
        with pytest.raises(ValueError):
            db.get_ephemeris_exact('2P', '500', epochs, source='horizons')

    def test_get_ephemeris_interp(self, db):
        jdc = 2458120.0
        jda, jdb = 2458119.5, 2458120.5
        eph = db.get_ephemeris(2, jda, jdb)
        a = RADec(eph[0]['ra'], eph[0]['dec'], unit='rad')
        b = RADec(eph[1]['ra'], eph[1]['dec'], unit='rad')
        test = util.spherical_interpolation(a, b, jda, jdb, jdc)

        eph, v = db.get_ephemeris_interp(2, [jdc])
        assert np.isclose(eph.separation(test), 0)

    def test_get_ephemeris_segments(self, db):
        jda, jdb = 2458119.5, 2458121.5
        db.add_ephemeris(1, '500', jda, jdb, step='1d', source='mpc',
                         cache=True)
        segments = db.get_ephemeris_segments()
        assert len(list(segments)) == 6

        segments = db.get_ephemeris_segments(start=jda, stop=jdb - 1)
        assert len(list(segments)) == 4

        segments = db.get_ephemeris_segments(objid=1, start=jda, stop=jdb - 1)
        assert len(list(segments)) == 2

    def test_get_found_by_id(self, db):
        c = db.get_found_by_id([1, 3], columns='*')
        assert len(c) == 2

        c = db.get_found_by_id([1, 3], columns='count()', generator=True)
        assert len(list(c)) == 2

        assert len(db.get_found_by_id([100])) == 0

    def test_get_found_date(self, db):
        start = 2458119.5
        stop = start + 60 / 86400

        c = db.get_found(start=start, stop=stop, columns='count()')[0]
        assert c[0] == 1

        c = db.get_found(start=start, stop=stop, columns='count()',
                         generator=True)
        assert next(c)[0] == 1

    def test_get_found_by_obsid(self, db):
        c = db.get_found_by_obsid([1])
        assert len(c) == 1

        c = list(db.get_found_by_obsid([2], generator=True))
        assert len(c) == 1

        c = db.get_found_by_obsid([5])
        assert len(c) == 0

    def test_get_found_by_object(self, db):
        c = db.get_found(obj=1)
        assert len(c) == 0

        c = list(db.get_found(obj=2, generator=True))
        assert len(c) == 3

    def test_get_objects(self, db):
        objid, desg = db.get_objects()
        assert len(objid) == 2
        assert 'C/1995 O1' in desg
        assert '2P' in desg

    def test_get_observations_by_date(self, db):
        obsids = db.get_observations_by_date(2458119.5, 2458121.5,
                                             columns='obsid')
        assert len(obsids) == N_tiles**2

    def test_get_observations_by_id(self, db):
        obsids = db.get_observations_by_id([1, 2, 3], generator=True)
        assert len(list(obsids)) == 3

        assert len(db.get_observations_by_id([100])) == 0

    def test_get_observations_overlapping(self, db):
        eph = db.get_ephemeris(2, None, None)

        ra = [eph[i]['ra'] for i in range(len(eph))]
        dec = [eph[i]['dec'] for i in range(len(eph))]

        epochs = [eph[i]['jd'] for i in range(len(eph))]
        start = min(epochs)
        stop = max(epochs)

        obs = db.get_observations_overlapping(
            ra=ra, dec=dec, start=start, stop=stop)

        # for N_tiles == 10, ephemeris will be in just one box
        assert len(obs) == 1

        obs = db.get_observations_overlapping(ra=ra, start=start, stop=stop)
        assert len(obs) == 4

        obs = db.get_observations_overlapping(dec=dec, start=start, stop=stop)
        assert len(obs) == 10

    def test_get_observations_overlapping_error(self, db):
        with pytest.raises(ValueError):
            db.get_observations_overlapping(box={}, ra=[1])

        with pytest.raises(ValueError):
            db.get_observations_overlapping(ra=[1])

        with pytest.raises(ValueError):
            db.get_observations_overlapping(dec=[1])

        with pytest.raises(ValueError):
            db.get_observations_overlapping(ra=[1, 2, 3], dec=[4, 5, 6, 7])

    def test_get_orbit_exact(self, db):
        epochs = [2458119.5]
        orb = db.get_orbit_exact(2, [2458119.5], cache=True)
        assert len(orb) == 1

    def test_resolve_objects(self, db):
        objid, desg = list(zip(*db.resolve_objects([1, '2P'])))
        assert objid[0] == 1
        assert objid[1] == 2
        assert desg[0] == 'C/1995 O1'
        assert desg[1] == '2P'

    def test_resolve_object(self, db):
        objid, desg = db.resolve_object(1)
        assert objid == 1
        assert desg == 'C/1995 O1'

        objid, desg = db.resolve_object('2P')
        assert objid == 2
        assert desg == '2P'

        objid, desg = db.resolve_object('1P')
        assert objid is None

    def test_resolve_object_fail(self, db):
        with pytest.raises(BadObjectID):
            db.resolve_object(23)
