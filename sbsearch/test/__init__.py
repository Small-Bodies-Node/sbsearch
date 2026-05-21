import pytest
import testing.postgresql
import sqlalchemy as sa
from ..sbsearch import SBSearch
from ..sbsdb import SBSDatabase

Postgresql = testing.postgresql.PostgresqlFactory(cache_initialized_db=True)


def get_url(postgresql):
    # use psycopg dialect, which will not be the default for sqlalchemy until
    # 2.1
    return postgresql.url().replace("postgresql://", "postgresql+psycopg://")


@pytest.fixture(name="sbs")
def fixture_sbs() -> SBSearch:
    with Postgresql() as postgresql:
        engine: sa.engine.Engine = sa.create_engine(get_url(postgresql))
        sessionmaker: sa.orm.sessionmaker = sa.orm.sessionmaker(bind=engine)
        with SBSearch(
            sessionmaker(), min_edge_length=0.01, max_edge_length=0.17
        ) as sbs:
            yield sbs


@pytest.fixture(name="db")
def fixture_db():
    with Postgresql() as postgresql:
        db = SBSDatabase.test_db(get_url(postgresql))
        yield db
        db.close()
