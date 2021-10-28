# Licensed with the 3-clause BSD license.  See LICENSE for details.

import pytest

import sqlalchemy as sa
from sqlalchemy.exc import ProgrammingError

from . import fixture_db, Postgresql
from ..sbsdb import SBSDatabase


class TestSBSDatabase:
    def test_init_session(self):
        with Postgresql() as postgresql:
            engine = sa.create_engine(postgresql.url())
            sessionmaker = sa.orm.sessionmaker(bind=engine)
            session = sessionmaker()
            db = SBSDatabase(session)
            db.create()
            db.session.execute('SELECT * FROM designation').fetchall()

    def test_verify(self, db):
        db.verify()

    def test_verify_missing_table(self, db):
        db.session.execute('DROP TABLE designation')
        with pytest.raises(ProgrammingError):
            db.session.execute('SELECT * FROM designation').fetchall()
        db.verify()
        db.session.execute('SELECT * FROM designation').fetchall()
