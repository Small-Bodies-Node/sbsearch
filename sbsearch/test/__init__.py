import pytest
import testing.postgresql
import sqlalchemy as sa
from ..sbsearch import SBSearch
from ..sbsdb import SBSDatabase

# Use `handler()` on initialize database


def handler(postgresql):
    engine = sa.create_engine(postgresql.url())
    with engine.connect() as connection:
        connection.execute('CREATE EXTENSION IF NOT EXISTS pg_trgm')


Postgresql = testing.postgresql.PostgresqlFactory(cache_initialized_db=True,
                                                  on_initialized=handler)


@pytest.fixture(name='sbs')
def fixture_sbs() -> SBSearch:
    with Postgresql() as postgresql:
        engine: sa.engine.Engine = sa.create_engine(postgresql.url())
        sessionmaker: sa.orm.sessionmaker = sa.orm.sessionmaker(bind=engine)
        with SBSearch(sessionmaker()) as sbs:
            yield sbs


@pytest.fixture(name='db')
def fixture_db():
    with Postgresql() as postgresql:
        db = SBSDatabase.test_db(postgresql.url())
        yield db
        db.close()
