# Licensed with the 3-clause BSD license.  See LICENSE for details.
import pytest
from logging import Logger
from itertools import repeat

import numpy as np
import astropy.units as u
import sqlalchemy as sa
from sbpy.data import Orbit

from .. import util, schema
from ..schema import Obj, Obs, Eph, Found
from ..util import RADec, FieldOfView
from ..db import SBDB
from ..exceptions import BadObjectID, NoEphemerisError, SourceNotFoundError
from .skytiles import N_tiles


@pytest.fixture
def db():
    db = SBDB.create_test_db()
    yield db


class TestSBDB:
    def test_verify_database(self, db):
        logger = Logger('test')
        db.verify_database(logger)

        metadata = sa.MetaData(bind=db.engine)
        metadata.reflect()
        for table in db.DB_NAMES:
            assert table in metadata.tables.keys()

        db.engine.execute('DROP TABLE eph')
        metadata.clear()
        metadata.reflect()
        assert 'eph' not in metadata.tables.keys()
        db.verify_database(logger)
        metadata.clear()
        metadata.reflect()
        assert 'eph' in metadata.tables.keys()

        class Test(schema.Base):
            __tablename__ = 'test'
            id = sa.Column(sa.Integer, primary_key=True)

        assert 'test' not in metadata.tables.keys()
        db.verify_database(logger, names=['test'])
        metadata.clear()
        metadata.reflect()
        assert 'test' in metadata.tables.keys()

    def test_create_test_db(self):
        # just exercise the code
        db = SBDB.create_test_db()
        assert db is not None

    def test_add_ephemeris_fixed(self, db):
        # already added via test database
        c = db.session.query(Eph).count()
        assert c == 3  # default number in test database

    def test_add_ephemeris_adaptable(self, db):
        c0 = db.session.query(Eph).count()
        db.add_ephemeris(2, '500', 2457799.5, 2457809.5, step=None,
                         source='mpc', cache=True)
        c1 = db.session.query(Eph).count()
        assert c1 - c0 == 36

    def test_add_found(self, db):
        rows = db.get_found().all()
        assert len(rows) == 3

        foundids = db.add_found_by_id(2, [3], '500', cache=True,
                                      update=False)
        assert len(foundids) == 0

        foundids = db.add_found_by_id(2, [3], '500', cache=True,
                                      update=True)
        assert len(foundids) == 1

    def test_add_found_by_id(self, db):
        # check test db found items
        rows = db.get_found().all()
        assert len(rows) == 3

        foundids = db.add_found_by_id(2, [3], '500', cache=True,
                                      update=False)
        assert len(foundids) == 0

        foundids = db.add_found_by_id(2, [3], '500', cache=True, update=True)
        assert len(foundids) == 1

    def test_add_object(self, db):
        # check test db object
        obj = db.get_objects().filter_by(desg="C/1995 O1").one()
        assert obj.objid == 1
        assert obj.desg == 'C/1995 O1'

    def test_add_object_error(self, db):
        with pytest.raises(ValueError):
            db.add_object(1)

    def test_add_observations(self, db):
        c = db.session.query(Obs).count()
        assert c == N_tiles**2

        # add rows
        obs = [Obs(
            jd_start=2450005.0,
            jd_stop=2450010.0,
            fov='SRID=40001;POLYGON((0 0, 1 1, 1 -1, -1 -1, -1 1, 0 0))',
        )]
        db.add_observations(obs)
        c = db.session.query(Obs).count()
        assert c == N_tiles**2 + 1

    def test_add_observations_with_source(self, db):
        class Survey(Obs):
            __tablename__ = 'survey'
            id = sa.Column(sa.Integer, primary_key=True)
            obsid = sa.Column(sa.Integer, sa.ForeignKey('obs.obsid'))

        with db.engine.connect() as con:
            schema.Base.metadata.create_all(con)

        fov = FieldOfView(RADec(np.random.rand(4, 2), unit='deg'))
        obs = [Survey(
            jd_start=2458300.5,
            jd_stop=2458300.51,
            fov=str(fov)
        )]
        db.add_observations(obs)

        c = db.session.query(Obs).count()
        assert c == N_tiles**2 + 1

        c = db.session.query(Survey).count()
        assert c == 1

    def test_clean_ephemeris(self, db):
        jda, jdb = 2458119.5, 2458121.5
        c = db.get_ephemeris(2, jda, jdb).count()
        assert c == 3
        c = db.clean_ephemeris(2, jda, jdb)
        assert c == 3
        c = db.get_ephemeris(2, jda, jdb).count()
        assert c == 0

    def test_clean_found(self, db):
        count = db.clean_found(2, None, None)
        assert count == 3

    # tested up to here
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
        c = db.get_found_by_id([1, 3]).count()
        assert c == 2

        assert db.get_found_by_id([100]).count() == 0

    def test_get_found_date(self, db):
        start = 2458119.5
        stop = start + 60 / 86400

        c = db.get_found(start=start, stop=stop).count()
        assert c == 1

    def test_get_found_by_obsid(self, db):
        c = db.get_found_by_obsid([1]).count()
        assert c == 0

    def test_get_found_by_object(self, db):
        c = db.get_found(obj=1).count()
        assert c == 0

        c = db.get_found(obj=2).count()
        assert c == 3

    def test_get_objects(self, db):
        objid, desg = db.get_objects()
        assert len(objid) == 2
        assert 'C/1995 O1' in desg
        assert '2P' in desg

    def test_get_observation_date_range(self, db):
        jd_range = db.get_observation_date_range()
        assert np.allclose(jd_range, [2458119.5, 2458119.5 + 50 / 1440])

        fov = str(FieldOfView(RADec(np.random.rand(4, 2))))
        obs = []
        obs.append(Obs(jd_start=2458300.5, jd_stop=2458300.51, fov=fov))
        obs.append(Obs(jd_start=2458100.5, jd_stop=2458100.51, fov=fov))

        jd_range = db.get_observation_date_range(source='test')
        assert np.allclose(jd_range, [2458119.5, 2458119.5 + 50 / 1440])

    def test_get_observation_date_range_error(self, db):
        with pytest.raises(SourceNotFoundError):
            jd_range = db.get_observation_date_range(source='blah')

    def test_get_observations_by_date(self, db):
        c = db.get_observations_by_date(
            2458119.5, 2458121.5, columns='obsid').count()
        assert c == N_tiles**2

    def test_get_observations_by_id(self, db):
        c = db.get_observations_by_id([]).count()
        assert c == 0

        c = db.get_observations_by_id([1, 2, 3]).count()
        assert c == 3

        c = db.get_observations_by_id([100]).count()
        assert c == 0

    def test_get_observations_containing(self, db):
        eph = db.get_ephemeris(2, None, None)
        coords = RADec.from_eph(eph)

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
