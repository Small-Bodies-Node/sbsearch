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
* [sbpy](https://github.com/NASA-Planetary-Science/sbpy) recent dev version

Optional packages:
* [pyoorb](https://github.com/oorb/oorb) for MPC Possible Comet Confirmation Page checking
* pytest for running the tests

## Testing
```
python setup.py build_ext --inplace
pytest sbsearch
```

## Contact

Maintained by [Michael S. P. Kelley](https://github.com/mkelley).  File an issue with your questions.

## References

Orbit integrations by OpenOrb: Granvik et al. Granvik, M., Virtanen, J., Oszkiewicz, D., Muinonen, K. (2009).  OpenOrb: Open-source asteroid orbit computation software including statistical ranging.  Meteoritics & Planetary Science 44(12), 1853-1861.
