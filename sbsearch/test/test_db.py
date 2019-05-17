# Licensed with the 3-clause BSD license.  See LICENSE for details.
import sqlite3
import pytest
from logging import Logger
from itertools import repeat

import numpy as np
import astropy.units as u
from sbpy.data import Orbit

from .. import util
from ..util import RADec
from ..db import SBDB
from ..exceptions import BadObjectID, NoEphemerisError, SourceNotFoundError
from .skytiles import sky_tiles, N_tiles


@pytest.fixture
def db():
    db = sqlite3.connect(':memory:', 5, 0, None, True, SBDB)
    db.verify_database(Logger('test'))
    db.add_object('C/1995 O1')
    objid = db.add_object('2P')
    db.add_alternate_desg(objid, 'Encke')

    obsids = range(N_tiles**2)
    start = 2458119.5 + np.arange(N_tiles**2) * 30 / 86400
    stop = start + 30 / 86400

    columns = [obsids, repeat('test'), start, stop, sky_tiles]
    db.add_observations(zip(*columns))
    db.add_ephemeris(2, '500', 2458119.5, 2458121.5, step='1d',
                     source='mpc', cache=True)
    db.add_found_by_id(2, [1, 2, 3], '500', cache=True)

    yield db
    db.close()


class TestSBDB:
    def test_verify_database(self, db):
        logger = Logger('test')
        db.verify_database(logger)
        c = db.execute("SELECT name FROM sqlite_master")
        existing_names = list([row[0] for row in c])
        count = sum([name in db.DB_NAMES for name in existing_names])
        assert count == len(db.DB_NAMES)

        db.execute('drop table eph')
        db.verify_database(logger)
        c = db.execute("SELECT name FROM sqlite_master")
        existing_names = list([row[0] for row in c])
        count = sum([name in db.DB_NAMES for name in existing_names])
        assert count == len(db.DB_NAMES)

        script = 'CREATE TABLE test(test);'
        db.verify_database(logger, names=['test'], script=script)
        c = db.execute("SELECT name FROM sqlite_master")
        existing_names = list([row[0] for row in c])
        assert 'test' in existing_names

    def test_make_test_db(self):
        # just exercise the code
        db = SBDB.make_test_db()
        assert db is not None

    def test_add_alternate_desg(self, db):
        objid = db.resolve_object('2P')[0]
        c = db.execute('SELECT crossid, objid FROM altobj WHERE desg=?',
                       ('Encke',)).fetchall()
        assert len(c) == 1
        assert c[0][0] == 1
        assert c[0][1] == objid

    def test_add_alternate_desg_error(self, db):
        objid = db.resolve_object('2P')[0]
        with pytest.raises(ValueError):
            db.add_alternate_desg(objid, '2P')

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

        observations = db.get_observations_by_id([3])
        foundids = db.add_found_by_id(2, observations, '500', cache=True,
                                      update=False)
        assert len(foundids) == 0

        foundids = db.add_found_by_id(2, observations, '500', cache=True,
                                      update=True)
        assert len(foundids) == 1

    def test_add_found_by_id(self, db):
        rows = db.execute('select * from found').fetchall()
        assert len(rows) == 3

        foundids = db.add_found_by_id(2, [3], '500', cache=True,
                                      update=False)
        assert len(foundids) == 0

        foundids = db.add_found_by_id(2, [3], '500', cache=True, update=True)
        assert len(foundids) == 1

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
        rows = [[100, 'test', 5, 10, [0, 0, 1, 1, 1, -1, -1, -1, -1, 1]]]
        db.add_observations(rows)
        c = db.execute('select count() from obs').fetchone()[0]
        assert c == N_tiles**2 + 1

    def test_add_observations_with_other_table(self, db):
        db.execute('''
        CREATE TABLE survey(obsid INTEGER PRIMARY KEY, a INTEGER)
        ''')

        points = np.random.rand(10)
        new_obs = [[None, 'survey', 2458300.5, 2458300.51, points]]
        other_cmd = 'INSERT INTO survey VALUES (last_insert_rowid(),?)'
        other_rows = [[5]]
        db.add_observations(new_obs, other_cmd=other_cmd,
                            other_rows=other_rows, logger=Logger('test'))

        c = db.execute('select count() from obs').fetchone()[0]
        assert c == N_tiles**2 + 1

        c = db.execute('select count() from survey').fetchone()[0]
        assert c == 1

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

    def test_get_alternates(self, db):
        encke = db.resolve_object('2P')[0]
        db.add_alternate_desg(encke, 'encke')

        hb = db.resolve_object('C/1995 O1')[0]
        db.add_alternate_desg(hb, 'Hale-Bopp')

        alternates = db.get_alternates()
        assert len(alternates[encke]) == 2
        assert 'Encke' in alternates[encke]
        assert 'encke' in alternates[encke]
        assert len(alternates[hb]) == 1
        assert 'Hale-Bopp' in alternates[hb]

    def test_get_ephemeris_date_range(self, db):
        jd_range = db.get_ephemeris_date_range()
        assert np.allclose(jd_range, [2458119.5, 2458121.5])

        jd_range = db.get_ephemeris_date_range(objids=[2])
        assert np.allclose(jd_range, [2458119.5, 2458121.5])

    def test_get_ephemeris_date_range_error(self, db):
        with pytest.raises(NoEphemerisError):
            jd_range = db.get_ephemeris_date_range([1])

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
        eph = db.get_ephemeris_exact('2P', '500', epochs, source='oorb',
                                     orbit=orbit)
        assert len(eph) == 3

    def test_get_ephemeris_exact_error(self, db):
        epochs = (2458119.5, 2458120.5, 2458121.5)
        with pytest.raises(ValueError):
            db.get_ephemeris_exact('2P', '500', epochs, source='horizons')

        epochs = (2458119.5, 2458119.5, 2458119.5)
        with pytest.raises(ValueError):
            db.get_ephemeris_exact('2P', '500', epochs, source='jpl')

    def test_get_ephemeris_exact_365(self, db):
        epochs = np.arange(365) + 2458119.5
        eph = db.get_ephemeris_exact('2P', '500', epochs, source='jpl',
                                     cache=True)
        assert len(eph) == 365

    def test_get_ephemeris_interp(self, db):
        jdc = 2458120.0
        jda, jdb = 2458119.5, 2458120.5
        eph = db.get_ephemeris(2, jda, jdb)
        a = RADec(eph[0][5], eph[0][6], unit='rad')
        b = RADec(eph[1][5], eph[1][6], unit='rad')
        test = util.spherical_interpolation(a, b, jda, jdb, jdc)

        eph, v = db.get_ephemeris_interp(2, [jdc])
        assert np.isclose(eph.separation(test), 0)

    def test_get_ephemeris_segments(self, db):
        jda, jdb = 2458119.5, 2458121.5
        db.add_ephemeris(1, '500', jda, jdb, step='1d', source='mpc',
                         cache=True)
        ephids, segments = db.get_ephemeris_segments()
        assert len(ephids) == 6
        for k, v in segments.items():
            assert len(v) == 6

        ephids, segments = db.get_ephemeris_segments(
            start=jda, stop=jdb - 1)
        assert len(ephids) == 4

        ephids, segments = db.get_ephemeris_segments(
            objid=1, start=jda, stop=jdb - 1)
        assert len(ephids) == 2

        ephids, segments = db.get_ephemeris_segments(
            objid=1, start=jda, stop=jdb - 1, vmax=1)
        assert len(ephids) == 0

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
        objid, desg, alternates = db.get_objects()
        assert len(objid) == 2
        assert 'C/1995 O1' in desg
        assert '2P' in desg
        assert 'Encke' in alternates

    def test_get_observation_date_range(self, db):
        jd_range = db.get_observation_date_range()
        assert np.allclose(jd_range, [2458119.5, 2458119.5 + 50 / 1440])

        points = np.random.rand(10)
        db.add_observations([[None, 'blah', 2458300.5, 2458300.51, points]])
        db.add_observations([[None, 'blah', 2458100.5, 2458100.51, points]])
        jd_range = db.get_observation_date_range(source='test')
        assert np.allclose(jd_range, [2458119.5, 2458119.5 + 50 / 1440])

    def test_get_observation_date_range_error(self, db):
        with pytest.raises(SourceNotFoundError):
            jd_range = db.get_observation_date_range(source='blah')

    def test_get_observations_by_date(self, db):
        obsids = db.get_observations_by_date(
            2458119.5, 2458121.5, columns='obsid')
        assert len(obsids) == N_tiles**2

    def test_get_observations_by_id(self, db):
        obs = db.get_observations_by_id([])
        assert len(obs) == 0

        obs = db.get_observations_by_id([1, 2, 3], generator=True)
        assert len(list(obs)) == 3

        assert len(db.get_observations_by_id([100])) == 0

    def test_get_observations_by_id_inner_join(self, db):
        db.execute('''
        CREATE TABLE survey(obsid INTEGER PRIMARY KEY, a INTEGER)
        ''')

        points = np.random.rand(10)
        new_obs = [[1001, 'survey', 2458300.5, 2458300.51, points]]
        other_cmd = 'INSERT INTO survey VALUES (last_insert_rowid(),?)'
        other_rows = [[5]]
        db.add_observations(new_obs, other_cmd=other_cmd,
                            other_rows=other_rows, logger=Logger('test'))

        obs = db.get_observations_by_id(
            [1001], inner_join=['survey USING (obsid)'])
        assert len(list(obs)) == 1

    def test_get_observations_near(self, db):
        eph = db.get_ephemeris(2, None, None)

        ra = [eph[i][5] for i in range(len(eph))]
        dec = [eph[i][6] for i in range(len(eph))]

        epochs = [eph[i][2] for i in range(len(eph))]
        start = min(epochs)
        stop = max(epochs)

        obs = db.get_observations_near(
            ra=ra, dec=dec, start=start, stop=stop)

        # for N_tiles == 10, ephemeris will be in just one box
        assert len(obs) == 1

        obs = db.get_observations_near(ra=ra, start=start, stop=stop)
        assert len(obs) == 4

        obs = db.get_observations_near(dec=dec, start=start, stop=stop,
                                       generator=True)
        assert len(list(obs)) == 10

    def test_get_observations_near_error(self, db):
        with pytest.raises(ValueError):
            db.get_observations_near(ra=[1])

        with pytest.raises(ValueError):
            db.get_observations_near(dec=[1])

        with pytest.raises(ValueError):
            db.get_observations_near(ra=[1, 2, 3], dec=[4, 5, 6, 7])

    def test_get_observations_near_box(self, db):
        # need a dummy table to test out inner_join
        db.execute('CREATE TABLE dummy(obsid INTEGER PRIMARY KEY, a FLOAT)')
        db.execute('INSERT INTO dummy SELECT obsid,1.0 FROM obs')

        ra = np.radians(np.linspace(0, 50, 30))
        dec = np.radians(np.linspace(0, 50, 30))
        x, y, z = util.rd2xyz(ra, dec)
        x = x.reshape(3, 10)
        y = y.reshape(3, 10)
        z = z.reshape(3, 10)

        box = {'x0': x.min(0), 'x1': x.max(0),
               'y0': y.min(0), 'y1': y.max(0),
               'z0': z.min(0), 'z1': z.max(0)}
        obs = db.get_observations_near_box(
            inner_join=['dummy USING (obsid)'], **box)
        assert len(obs) == 6

    def test_get_observations_near_box_error(self, db):
        with pytest.raises(ValueError):
            db.get_observations_near_box()

    def test_get_orbit_exact(self, db):
        orb = db.get_orbit_exact(2, [2458119.5], cache=True)
        assert len(orb) == 1

        epochs = {'start': '2018-01-01', 'stop': '2018-01-03', 'step': '1d'}
        orb = db.get_orbit_exact(2, epochs, cache=True)
        assert len(orb) == 3

    def test_get_orbit_exact_365(self, db):
        epochs = 2458119.5 + np.arange(365)
        orb = db.get_orbit_exact(2, epochs, cache=True)
        assert len(orb) == 365

        epochs = {'start': '2018-01-01', 'stop': '2018-01-03', 'step': '1d'}
        orb = db.get_orbit_exact(2, epochs, cache=True)
        assert len(orb) == 3

    def test_get_orbit_exact_error(self, db):
        epochs = [2458119.5] * 3
        with pytest.raises(ValueError):
            orb = db.get_orbit_exact(2, epochs, cache=True)

    def test_remove_alternate_designation(self, db):
        c = db.remove_alternate_designation('Encke')
        assert c == 1

        c = db.remove_alternate_designation('Encke')
        assert c == 0

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

        # test altobj
        objid, desg = db.resolve_object('Encke')
        assert objid == 2
        assert desg == '2P'

    def test_update_object(self, db):
        objid = db.resolve_object('2P')[0]
        db.update_object(objid, '2P/Encke')

        new_desg = db.resolve_object('2P/Encke')[1]
        assert new_desg == '2P/Encke'

        alternates = db.get_alternates(objid)
        assert len(alternates) == 2
        assert '2P' in alternates
        assert 'Encke' in alternates  # already defined

        # back to 2P
        db.update_object(objid, '2P')
        new_desg = db.resolve_object('2P/Encke')[1]
        assert new_desg == '2P'

        alternates = db.get_alternates(objid)
        assert len(alternates) == 2
        assert '2P/Encke' in alternates
        assert 'Encke' in alternates

    def test_update_object_by_objid_error(self, db):
        with pytest.raises(TypeError):
            db.update_object('2P', '2P/Encke')

    def test_resolve_object_fail(self, db):
        with pytest.raises(BadObjectID):
            db.resolve_object(23)
