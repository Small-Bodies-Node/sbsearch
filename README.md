# sbsearch v0.2.0-dev
Search for specific small Solar System bodies in astronomical surveys.

`sbsearch` is designed for efficient searching of large amounts of wide-field data.  The guiding principle is to execute a fast and approximate search to narrow down the list of images and objects needed for a more-precise search.   The search is based on ephemerides from the Minor Planet Center or JPL Horizons.  Ephemerides for objects commonly searched for can be stored and re-used.

## Status

2019 Feb 09: `sbsearch` is functional, but in development and therefore the API can still change.  Use the v0.1 branch for a more stable code.

## Requirements

* Python 3.5+
* requests
* cython
* sqlite 3+
* wget
* astropy 2.0+
* [astroquery](https://astroquery.readthedocs.io/en/latest/) 0.3.9+
* [sbpy](https://github.com/NASA-Planetary-Science/sbpy) 0.1

Optional packages:
* [pyoorb](https://github.com/oorb/oorb) for MPC Possible Comet Confirmation Page checking
* pytest for running the tests

## Usage

`sbsearch` can be used as is, but generally you'll want to add survey specific metadata.  A few columns are already defined: filter, seeing, airmass, and maglimit.  To add other metadata, subclass the `Obs` object for your survey, and define the necessary attributes:

``` python
from sqlalchemy import Integer, Float, String, ForeignKey
from sbsearch.schema import Obs
class ZTF(Obs):
    __tablename__ = 'ztf'
    pid = Column(Integer, primary_key=True)
    obsid = Column(ForeignKey("obs.obsid"))
	obsdate = Column(String(32))
	infobits = Column(Integer)
	field = Column(Integer)
	ccdid = Column(Integer)

    __mapper_args__ = {
        'polymorphic_identity': 'ztf',
    }
```

If this is defined before the `SBSearch` object is initialized, then your table will also be created and the new survey object may be used in place of `Obs` for inserting observations:

``` python
sbs = SBSearch()
ztf_obs = ZTF(
    jd_start=2458606.147528218,
	jd_stop=2458606.147630901,
	filter='zr',
	seeing=2.2,
	airmass=1.4,
	maglimit=20.3
	pid=12345,
	obsdate='2019-05-02',
	infobits=0,
	field=512,
	ccdid=11)
ztf_obs.coords_to_fov((0, 1), (1, 1), (1, 0), (0, 0), (0.5, 0.5))
sbs.add_observation(ztf_obs)
```

## Testing
```
python setup.py build_ext --inplace
pytest sbsearch
```

## Contact

Maintained by [Michael S. P. Kelley](https://github.com/mkelley).  File an issue with your questions.

## References

Orbit integrations by OpenOrb: Granvik et al. Granvik, M., Virtanen, J., Oszkiewicz, D., Muinonen, K. (2009).  OpenOrb: Open-source asteroid orbit computation software including statistical ranging.  Meteoritics & Planetary Science 44(12), 1853-1861.
