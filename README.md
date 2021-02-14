# sbsearch v2.0.0-dev
Search for specific small Solar System bodies in astronomical surveys.

`sbsearch` is designed for efficient searching of large amounts of wide-field data.  The guiding principle is to execute a fast and approximate search to narrow down the list of images and objects needed for a more-precise search.  The search is based on ephemerides from the Minor Planet Center or JPL Horizons.  Ephemerides for objects commonly searched for can be stored and re-used.

v2 is a complete re-write, primarily to replace PostGIS with s2geometry and enable searches considering ephemeris uncertainties.  The code is conceptually similar to but incompatible with previous versions.

## Requirements

* Python 3.6+
* [s2geometry](s2geometry.io)
* cython
* [SQLAlchemy](https://www.sqlalchemy.org/) 1.3
* A database backend, e.g., sqlite3 or PostgresSQL.  A database dialect for SQLAlchemy may also be needed, e.g., psycopg2 for PostgreSQL.
* astropy 3.3+
* [astroquery](https://astroquery.readthedocs.io/en/latest/) 0.4.1+
* [sbpy](https://github.com/NASA-Planetary-Science/sbpy) 0.2.2

Optional packages:
* pytest, pytest-cov for running the tests


## Testing

Build the Cython extensions in place, and run the tests.  For example:
```
python3 setup.py build_ext --inplace
pytest sbsearch --cov=sbsearch --cov-report=html
```

Or, to find libs2 in a virtual environement:
```
LDFLAGS="-L$VIRTUAL_ENV/lib -Wl,-rpath=$VIRTUAL_ENV/lib" python3 setup.py build_ext --inplace
pytest sbsearch --cov=sbsearch --cov-report=html
```


## Usage

### Survey-specific metadata

`sbsearch` can be used as is, but generally you'll want to add survey specific metadata.  A few columns are already defined: filter, seeing, airmass, and maglimit.  To add other metadata and survey specific parameters (name, observatory location), subclass the `Observation` object for your survey, and define the necessary attributes.  The object ``sbsearch.model.UnspecifiedSurvey`` can be used as an example.

```python
from sqlalchemy import Integer, Float, String, ForeignKey
from sbsearch.model import Observation
class ZTF(Observation):
    __tablename__ = 'ztf'
    __obscode__ = 'I41'  # ZTF's IAU observatory code
    pid = Column(Integer, primary_key=True)
    observation_id = Column(
        Integer, ForeignKey('observation.observation_id', onupdate='CASCADE',
                            ondelete='CASCADE'))
    infobits = Column(Integer)
    field = Column(Integer)
    ccdid = Column(Integer)
    qid = Column(Integer)

    __mapper_args__ = {
        'polymorphic_identity': 'ztf',
    }
```
With this object defined, the database will be updated the next time the `SBSearch` object is initialized.  The new survey object may be used in place of `Observation` for inserting observations:

``` python
sbs = SBSearch()
ztf_obs = ZTF(
    mjd_start=58605.647528218,
    mjd_stop=58605.647630901,
    filter='zr',
    seeing=2.2,
    airmass=1.4,
    maglimit=20.3
    pid=12345,
    obsdate='2019-05-02',
    infobits=0,
    field=512,
    ccdid=11,
    qid=1)
ztf_obs.set_fov((0, 1, 1, 0) (1, 1, 0, 0))
sbs.add_observation(ztf_obs)
```

## Contact

Maintained by [Michael S. P. Kelley](https://github.com/mkelley).  File an issue with your questions.

## References


## Developer notes
### s2geometry
gtest is supposed to be optional and there is a PR to fix that.  Until it is merged:
```
wget https://patch-diff.githubusercontent.com/raw/google/s2geometry/pull/78.patch
git apply --stat 78.patch
git apply --check 78.patch
git apply 78.patch
```

To build s2geometry to a virtual environment directory:
```
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV ..
make
make install
```
