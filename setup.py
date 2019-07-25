#!/usr/bin/env python
from setuptools import setup, find_packages


if __name__ == "__main__":
    setup(name='sbsearch',
          version='1.0.2',
          description='Find small Solar System bodies in astronomical surveys.',
          author="Michael S. P. Kelley",
          author_email="msk@astro.umd.edu",
          url="https://github.com/mkelley/sbsearch",
          packages=find_packages(),
          requires=['numpy', 'astropy', 'astroquery', 'sbpy',
                    'sqlalchemy', 'geoalchemy2'],
          setup_requires=['pytest-runner'],
          tests_require=['pytest'],
          license='BSD',
          )
