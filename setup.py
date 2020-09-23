#!/usr/bin/env python
from setuptools import setup, find_packages


if __name__ == "__main__":
    setup(name='sbsearch',
          version='1.1.3',
          description='Find small Solar System bodies in astronomical surveys.',
          author="Michael S. P. Kelley",
          author_email="msk@astro.umd.edu",
          url="https://github.com/mkelley/sbsearch",
          packages=find_packages(),
          install_requires=['numpy>=1.13', 'astropy<4.0', 'astroquery>=0.4.dev5744', 'sbpy>=0.2',
                            'sqlalchemy>=1.3', 'geoalchemy2'],
          setup_requires=['pytest-runner'],
          tests_require=['pytest'],
          license='BSD',
          )
