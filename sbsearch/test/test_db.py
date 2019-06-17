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
    db.close()


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

    def test_create_test_db(self, db):
        # just exercise the code
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
        logger = Logger('test')

        c = db.session.query(Obs).count()
        assert c == N_tiles**2

        # add rows
        obs = [Obs(
            jd_start=2450005.0,
            jd_stop=2450010.0,
            fov='SRID=40001;POLYGON((0 0, 1 1, 1 -1, -1 -1, -1 1, 0 0))',
        )]
        db.add_observations(obs, logger=logger)
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
            jd_range = db.get_ephemeris_date_range([3])

    def test_get_ephemeris(self, db):
        eph = db.get_ephemeris(2, 2458119.5, 2458121.5).all()
        assert len(eph) == 3

    def test_get_ephemeris_all(self, db):
        eph = db.get_ephemeris(None, None, None).all()
        assert len(list(eph)) == 3

    def test_get_ephemeris_interp(self, db):
        jdc = 2458120.0
        jda, jdb = 2458119.5, 2458120.5
        eph = db.get_ephemeris(2, jda, jdb).all()
        a = RADec(eph[0].ra, eph[0].dec, unit='rad')
        b = RADec(eph[1].ra, eph[1].dec, unit='rad')
        test = util.spherical_interpolation(a, b, jda, jdb, jdc)

        eph, v = db.get_ephemeris_interp(2, [jdc])
        assert np.isclose(eph.separation(test).value, 0)

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
        rows = db.get_objects().all()
        objid = [r.objid for r in rows]
        desg = [r.desg for r in rows]
        assert len(objid) == 2
        assert 'C/1995 O1' in desg
        assert '2P' in desg

    def test_get_observation_date_range(self, db):
        jd_range = db.get_observation_date_range()
        assert np.allclose(jd_range, [2458119.5, 2458119.5 + 50 / 1440])

        fov = str(FieldOfView(RADec(np.random.rand(4, 2), unit='deg')))
        obs = []
        obs.append(Obs(jd_start=2458300.5, jd_stop=2458300.51, fov=fov))
        obs.append(Obs(jd_start=2458100.5, jd_stop=2458100.51, fov=fov))

        jd_range = db.get_observation_date_range(source='unspecified')
        assert np.allclose(jd_range, [2458119.5, 2458119.5 + 50 / 1440])

    def test_get_observation_date_range_error(self, db):
        with pytest.raises(SourceNotFoundError):
            jd_range = db.get_observation_date_range(source='blah')

        obs = db.get_observations_by_date().delete()
        with pytest.raises(SourceNotFoundError):
            jd_range = db.get_observation_date_range()

    def test_get_observations_by_id(self, db):
        c = db.get_observations_by_id([]).count()
        assert c == 0

        c = db.get_observations_by_id([1, 2, 3]).count()
        assert c == 3

        c = db.get_observations_by_id([100]).count()
        assert c == 0

    def test_get_observations_by_date(self, db):
        c = db.get_observations_by_date(2458119.5, 2458121.5).count()
        assert c == N_tiles**2

    def test_get_observations_covering(self, db):
        eph = db.get_ephemeris(2, None, None).all()
        coords = RADec.from_eph(eph)
        fov = str(FieldOfView(RADec.from_eph(eph)))

        epochs = [e.jd for e in eph]
        start = min(epochs)
        stop = max(epochs)

        nobs = (db.get_observations_covering(fov, start=start, stop=stop)
                .count())

        # for N_tiles == 10, ephemeris will fully be in just one box
        assert nobs == 1

    def test_get_observations_intersecting(self, db):
        eph = db.get_ephemeris(2, None, None).all()
        coords = RADec.from_eph(eph)
        fov = str(FieldOfView(RADec.from_eph(eph)))

        epochs = [e.jd for e in eph]
        start = min(epochs)
        stop = max(epochs)

        nobs = (db.get_observations_intersecting(fov, start=start, stop=stop)
                .count())

        # for N_tiles == 10, ephemeris will fully be in just one box
        assert nobs == 1

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
