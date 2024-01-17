# sbsearch v0.1.12

Search for specific small Solar System bodies in astronomical surveys.

`sbsearch` is designed for efficient searching of large amounts of wide-field data.  The guiding principle is to execute a fast and approximate search to narrow down the list of images and objects needed for a more-precise search.   The search is based on ephemerides from the Minor Planet Center or JPL Horizons.  Ephemerides for objects commonly searched for can be stored and re-used.

## Status

2018 Feb 09: v0.1 has been used with [ZChecker](https://github.com/mkelley/zchecker) since Nov 2018 with good results.

## Requirements

* Python 3.8+
* requests
* cython
* sqlite
* wget
* astropy 4
* numpy >=1.17<1.23
* [astroquery](https://astroquery.readthedocs.io/) 0.4.6
* [sbpy](https://sbpy.readthedocs.io/) 0.3

Optional packages:

* [pyoorb](https://github.com/oorb/oorb) for MPC Possible Comet Confirmation Page checking
* pytest for running the tests

## Acknowledgments

Support for SBSearch was provided by the NASA/University of Maryland/Minor Planet Center Augmentation through the NASA Planetary Data System Cooperative Agreement NNX16AB16A.

## Testing

`python setup.py build_ext --inplace`

`pytest sbsearch`

## Contact

File an issue or contact Mike Kelley ([@mkelley](https://github.com/mkelley)) with questions.
