#!/usr/bin/env python
from setuptools import setup, find_packages


if __name__ == "__main__":
    setup(name='sbsearch',
          version='1.1.6',
          description='Find small Solar System bodies in astronomical surveys.',
          author="Michael S. P. Kelley",
          author_email="msk@astro.umd.edu",
          url="https://github.com/mkelley/sbsearch",
          packages=find_packages(),
          install_requires=['numpy>=1.17', 'astropy>4<5', 'astroquery>=0.4.4.dev7007', 'sbpy>0.2.2',
                            'sqlalchemy>=1.3<1.4', 'geoalchemy2'],
          setup_requires=['pytest-runner'],
          tests_require=['pytest'],
          license='BSD',
          )
