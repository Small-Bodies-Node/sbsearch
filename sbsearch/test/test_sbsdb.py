# Licensed with the 3-clause BSD license.  See LICENSE for details.

import pytest

import sqlalchemy as sa
from sqlalchemy.exc import OperationalError

from ..sbsdb import SBSDatabase


@pytest.fixture
def db():
    db = SBSDatabase.test_db()
    yield db
    db.close()


class TestSBSDatabase:
    def test_init_session(self):
        engine = sa.create_engine('sqlite://')
        sessionmaker = sa.orm.sessionmaker(bind=engine)
        session = sessionmaker()
        db = SBSDatabase(session)
        db.create()
        db.session.execute('SELECT * FROM designation').fetchall()

    def test_verify(self, db):
        db.verify()

    def test_verify_missing_table(self, db):
        db.session.execute('DROP TABLE designation')
        with pytest.raises(OperationalError):
            db.session.execute('SELECT * FROM designation').fetchall()
        db.verify()
        db.session.execute('SELECT * FROM designation').fetchall()
