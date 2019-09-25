#!/usr/bin/env python
from setuptools import setup, find_packages
from Cython.Build import cythonize


if __name__ == "__main__":
    setup(name='sbsearch',
          version='0.1.6',
          description='Find small Solar System bodies in astronomical surveys.',
          author="Michael S. P. Kelley",
          author_email="msk@astro.umd.edu",
          url="https://github.com/mkelley/sbsearch",
          packages=find_packages(),
          requires=['numpy', 'astropy', 'sbpy', 'requests'],
          setup_requires=['pytest-runner'],
          tests_require=['pytest'],
          ext_modules=cythonize('sbsearch/interior.pyx'),
          license='BSD',
          )
