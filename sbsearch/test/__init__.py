import pytest
import testing.postgresql
import sqlalchemy as sa
from ..sbsearch import SBSearch
from ..sbsdb import SBSDatabase

Postgresql = testing.postgresql.PostgresqlFactory(cache_initialized_db=True)


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
