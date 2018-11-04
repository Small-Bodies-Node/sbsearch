# Licensed with the 3-clause BSD license.  See LICENSE for details.
import sqlite3
from logging import Logger

import pytest

from ..db import SBDB
from ..config import Config


@pytest.fixture
def db():
    db = sqlite3.connect(':memory:', 5, 0, None, True, SBDB)
    db.verify_tables(Logger('test'))
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
        objid = db.add_object('2P')
        db.add_ephemeris(objid, '500', 2458119.5, 2458121.5, step='1d',
                         source='mpc', cache=True)
        c = db.execute('select count() from eph').fetchone()[0]
        assert c == 3

    def test_add_ephemeris_mpc_variable(self, db):
        objid = db.add_object('2P')
        db.add_ephemeris(objid, '500', 2457799.5, 2457809.5, step=None,
                         source='mpc', cache=True)
        c = db.execute('select count() from eph').fetchone()[0]
        assert c == 36

    def test_add_object(self, db):
        objid = db.add_object('C/1995 O1')
        row = db.execute('select * from obj').fetchone()
        assert row[0] == objid
        assert row[1] == 'C/1995 O1'
