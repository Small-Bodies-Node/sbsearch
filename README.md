# sbsearch

[![CI Tests](https://github.com/Small-Bodies-Node/sbsearch/actions/workflows/ci-tests.yml/badge.svg)](https://github.com/Small-Bodies-Node/sbsearch/actions/workflows/ci-tests.yml)
[![codecov](https://codecov.io/gh/Small-Bodies-Node/sbsearch/graph/badge.svg?token=OOD89OKYA2)](https://codecov.io/gh/Small-Bodies-Node/sbsearch)

Search for specific small Solar System bodies in astronomical surveys.

`sbsearch` is designed for efficient searching of large amounts of wide-field data. The guiding principle is to execute a fast and approximate search to narrow down the list of images and objects needed for a more-precise search. The search is based on ephemerides from the Minor Planet Center or JPL Horizons. Ephemerides for objects commonly searched for can be stored and re-used.

v2 is a complete re-write, replacing PostGIS with the S2 library, and enabling areal searches, e.g., ephemerides with uncertainties. The code is conceptually similar to but incompatible with previous versions.

## Requirements

- Python 3.8+
- [s2geometry](http://s2geometry.io) v0.10.0
- Cython
- [SQLAlchemy](https://www.sqlalchemy.org/) 1.3
- PostgresSQL. A database dialect for SQLAlchemy may also be needed, e.g., psycopg.
- astropy >=4.3
- [astroquery](https://astroquery.readthedocs.io/en/latest/) 0.4.4dev7007+
- [sbpy](https://github.com/NASA-Planetary-Science/sbpy) >0.3.0

Optional packages:

- pytest, coverage, testing.postgresql and submodules for running the tests

## Usage

### Survey-specific metadata

`sbsearch` is intended to be used as a software dependency. It is up to the user to add survey specific metadata. A few columns are already defined, e.g., mjd_start, mjd_stop, filter, seeing, airmass, and maglimit. See the `Observation` class in `sbsearch.model.core.py` for all attributes. To add other metadata and survey specific parameters (name, observatory location), subclass the `Observation` object for your survey, and define the necessary attributes. The file `sbsearch.model.example_survey.py` should be used as an example. The [`catch`](https://github.com/Small-Bodies-Node/catch) program may also be referenced as an example.

```python
class ZTF(Observation):
    __tablename__ = 'ztf'
    __data_source_name__ = 'Zwicky Transient Facility'
    __obscode__ = 'I41'  # MPC observatory code

    source_id = Column(BigInteger, primary_key=True)
    observation_id = Column(BigInteger,
                            ForeignKey('observation.observation_id',
                                       onupdate='CASCADE',
                                       ondelete='CASCADE'),
                            nullable=False,
                            index=True)

    product_id = Column(String(64), doc='Archive product id',
                        unique=True, index=True)
    infobits = Column(Integer)
    field = Column(Integer)
    ccdid = Column(Integer)
    qid = Column(Integer)

    __mapper_args__ = {
        'polymorphic_identity': 'ztf'
    }
```

With this object defined, the database will be updated the next time the `SBSearch` object is initialized. The new survey object is used to insert observations:

```python
sbs = SBSearch()
obs = ZTF(
    mjd_start=58605.647528218,
    mjd_stop=58605.647630901,
    filter='zr',
    seeing=2.2,
    airmass=1.4,
    maglimit=20.3
    product_id=12345,
    infobits=0,
    field=512,
    ccdid=11,
    qid=1)
# set the observation's foot print on the sky:
ztf_obs.set_fov([0, 1, 1, 0], [1, 1, 0, 0])
sbs.add_observation(obs)
```

### Database maintenance

- After deleting any observations, the spatial index must be updated with

  ```sql
  REINDEX INDEX ix_observation_spatial_terms
  ```

- `VACUUM ANALYZE` may be useful after reindexing or adding new observations.

## Testing

Install testing dependencies, either manually or by pip installing the package
with `[test]`, e.g.,:

```bash
pip install -e .[test]
```

If installing dependencies manually, then build the Cython extensions in place:

```bash
python3 setup.py build_ext --inplace
```

Run the tests. For example:

```bash
pytest sbsearch
```

If libs2 is installed in your virtual environment, you may need:

```bash
LDFLAGS="-L$VIRTUAL_ENV/lib -Wl,-rpath=$VIRTUAL_ENV/lib" \
  CXXFLAGS="-I$VIRTUAL_ENV/include" \
  python3 setup.py build_ext --inplace
pytest sbsearch
```

Tests that require remote data (i.e., ephemerides) are skipped by default. To
run those tests:

```bash
pytest sbsearch --remote-data
```

Check test code coverage after running pytest by browsing `htmlcov/index.html`.

## Contact

Maintained by [Michael S. P. Kelley](https://github.com/mkelley). File an issue with your questions.

## Developer notes

### s2geometry

Requires: python3.x-dev, swig

s2 requires Abseil compiled with PIC enabled and the supporting the same C++ version. The following is an example the will compile Abseil and s2 to a Python virtual environment:

```bash
S2_TAG=v0.10.0
wget https://github.com/google/s2geometry/archive/refs/tags/${S2_TAG}.tar.gz
tar xzf ${S2_TAG}.tar.gz

ABSEIL_TAG=20220623.1
wget https://github.com/abseil/abseil-cpp/archive/refs/tags/${ABSEIL_TAG}.tar.gz
tar xzf ${ABSEIL_TAG}.tar.gz

cd abseil-cpp-${ABSEIL_TAG}
mkdir -p build
cd build
cmake -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_CXX_STANDARD=11 -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV -DABSL_ENABLE_INSTALL=ON -DABSL_PROPAGATE_CXX_STD=ON ..
make -j $(nproc)
make install
cd ../..
export LDFLAGS="-L$VIRTUAL_ENV/lib -Wl,-rpath=$VIRTUAL_ENV/lib"
export CXXFLAGS="-I$VIRTUAL_ENV/include"
export LD_LIBRARY_PATH=$VIRTUAL_ENV/lib

cd s2geometry-${S2_TAG:1}
mkdir -p build
cd build
cmake -DWITH_PYTHON=ON -DCMAKE_PREFIX_PATH=$VIRTUAL_ENV -DCMAKE_CXX_STANDARD=11 -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV -Wno-dev ..
make -j $(nproc)
make install
```

### Simulated data set

A large database based on a simulated query may be created and tested with the `test-extras/big-query.py` script. This requires a PostgreSQL database named `big_query_test` by default (`createdb big_query_test`).
