#!/usr/bin/env python
from setuptools import Extension, setup, find_packages
from Cython.Build import cythonize

if __name__ == "__main__":
    setup(
        name='sbsearch',
        version='2.0.0-dev',
        description='Find small Solar System bodies in astronomical surveys.',
        author="Michael S. P. Kelley",
        author_email="msk@astro.umd.edu",
        license='BSD',
        url="https://github.com/Small-Bodies-Node/sbsearch",
        packages=find_packages(),
        install_requires=['astropy>=3.2,<4', 'astroquery>=0.4.1', 'sbpy>=0.2.2',
                          'numpy<1.20', 'sqlalchemy>=1.3'],
        setup_requires=['pytest-runner', 'cython'],
        tests_require=['pytest'],
        ext_modules=cythonize(
            Extension("sbsearch.spatial", ["sbsearch/spatial.pyx"],
                      language="c++", libraries=['s2'])
        )
    )
