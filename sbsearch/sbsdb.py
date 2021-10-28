# Licensed with the 3-clause BSD license.  See LICENSE for details.
import logging
from typing import Type, TypeVar, Union

import sqlalchemy as sa
from sqlalchemy.orm import Session
from sqlalchemy.engine import Engine, reflection

from . import model

__all__ = ['SBSDatabase']


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
    db = SBSDatabase('postgresql://@/catch')

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
            self.engine = sa.create_engine(url_or_session, *args)
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

        """

        metadata: sa.MetaData = sa.MetaData()
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

        self.session.commit()

    def create_spatial_index(self):
        """Create the spatial term index.

        Generally VACUUM ANALZYE after this.

        """
        self.session.execute('''
        CREATE INDEX IF NOT EXISTS ix_observation_spatial_terms
        ON observation
        USING GIN (spatial_terms gin_trgm_ops);
        ''')
        self.session.commit()

    def drop_spatial_index(self):
        """Drop the spatial term index.

        Use this before inserting many observations.

        """
        self.session.execute('''
        DROP INDEX IF EXISTS ix_observation_spatial_terms;
        ''')
        self.session.commit()

    def create(self):
        model.Base.metadata.create_all(self.engine)
        self.create_spatial_index()

    @classmethod
    def test_db(cls: Type[SBSD], url: str) -> SBSD:
        from .target import MovingTarget

        db: SBSD = cls(url)
        db.create()

        MovingTarget('1P', db).add()
        target: MovingTarget = MovingTarget(
            'C/1995 O1', db,
            secondary_designations=['Hale-Bopp']
        ).add()

        return db
