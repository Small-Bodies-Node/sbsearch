# Licensed with the 3-clause BSD license.  See LICENSE for details.

try:
    from .version import version as __version__
except ImportError:
    __version__ = ''

from .sbsearch import *  # noqa

