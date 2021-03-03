# Licensed with the 3-clause BSD license.  See LICENSE for details.
import warnings
import logging
from typing import Type, TypeVar, Union

import sqlalchemy as sa
from sqlalchemy.orm import Session
from sqlalchemy.engine import Engine
from sqlalchemy.exc import ProgrammingError, SAWarning
import sqlite3
try:
    from psycopg2.extensions import connection as psycopg2_connection
    import psycopg2.extras
except ImportError:
    class psycopg2_connection:
        pass

from . import model

__all__ = ['SBSDatabase']


@sa.event.listens_for(Engine, "connect")
def set_sqlite_pragma(dbapi_connection, connection_record) -> None:
    if isinstance(dbapi_connection, sqlite3.Connection):
        cursor = dbapi_connection.cursor()
        cursor.execute("PRAGMA foreign_keys=ON")
        cursor.close()
    elif isinstance(dbapi_connection, psycopg2_connection):
        cursor = dbapi_connection.cursor()
        psycopg2.extras.register_range('floatrange', model.FloatRange, cursor,
                                       globally=True)
        cursor.close()


SBSD = TypeVar('SBSD', bound='SBSDatabase')


class SBSDatabase:
    """Database interface for SBSearch.


    Parameters
    ----------
    url_or_session : string or sqlalchemy Session
        The sqlalchemy-formatted database URL or a sqlalchemy session
        to use.

    *args
        `sqlalchemy.create_engine` arguments.


    Examples
    --------
    db = SBSDatabase('sqlite://')

    """

    def __init__(self, url_or_session: Union[str, Session], *args,
                 logger_name: str = 'SBSearch'):
        self.session: Session
        self.sessionmaker: Union[Session, None]
        self.engine: Engine
        if isinstance(url_or_session, Session):
            self.session = url_or_session
            self.engine = self.session.get_bind()
            self.sessionmaker = None
        else:
            self.engine = sa.create_engine(url_or_session, echo=False, *args)
            self.sessionmaker = sa.orm.sessionmaker(bind=self.engine)
            self.session = self.sessionmaker()

        self.logger: logging.Logger = logging.getLogger(logger_name)

    def __del__(self):
        self.close()

    def close(self):
        if self.session.is_active:
            self.session.close()

    def verify(self):
        """Verify SBSearch tables.

        Note: metadata.reflect will raise an exception when a table is
        missing, but another existing table refers to it (e.g., via
        foreign key).  To fix, drop the other table (preferred), or 
        manually create the missing table.

        Only indexes defined in sbsearch.model.INDICES are checked.

        """

        metadata: sa.MetaData = sa.MetaData()
        with warnings.catch_warnings():
            # suppress: SAWarning: Skipped unsupported reflection of
            # expression-based index ix_observation_mjd_index
            warnings.simplefilter("ignore", SAWarning)
            metadata.reflect(self.engine)

        missing: bool = False
        name: str
        for name in model.Base.metadata.tables.keys():
            if name not in metadata.tables.keys():
                missing = True
                self.logger.error('{} is missing from database'.format(name))

        if missing:
            self.create()
            self.logger.info('Created database tables.')

        for name in model.INDICES:
            try:
                getattr(model, name).create(self.engine)
            except ProgrammingError:
                pass

    def create(self):
        model.Base.metadata.create_all(self.engine)

    @classmethod
    def test_db(cls: Type[SBSD]) -> SBSD:
        from .target import MovingTarget

        db: SBSD = cls('sqlite://')
        db.create()

        MovingTarget('1P', db).add()
        target: MovingTarget = MovingTarget(
            'C/1995 O1', db,
            secondary_designations=['Hale-Bopp']
        ).add()

        return db
