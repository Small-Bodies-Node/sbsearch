from setuptools import Extension

def get_extensions():
    return [Extension("sbsearch.spatial", ["sbsearch/spatial.pyx"],
                      language="c++", libraries=['s2'])]
