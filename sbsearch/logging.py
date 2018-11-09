# Licensed under a 3-clause BSD style license - see LICENSE.rst
import sys
import logging
import numpy as np
from astropy.time import Time


def setup(filename='sbsearch.log', name='SBSearch'):
    logger = logging.Logger(name)
    logger.setLevel(logging.DEBUG)

    # This test allows logging to work when it is run multiple times
    # from ipython
    if len(logger.handlers) == 0:
        formatter = logging.Formatter('%(levelname)s: %(message)s')

        console = logging.StreamHandler(sys.stdout)
        console.setLevel(logging.DEBUG)
        console.setFormatter(formatter)
        logger.addHandler(console)

        logfile = logging.FileHandler(filename)
        logfile.setLevel(logging.INFO)
        logfile.setFormatter(formatter)
        logger.addHandler(logfile)

    logger.info('#' * 70)
    logger.info(Time.now().iso + 'Z')
    logger.info('Command line: ' + ' '.join(sys.argv[1:]))
    for handler in logger.handlers:
        if hasattr(handler, 'baseFilename'):
            logger.info('Logging to ' + handler.baseFilename)

    return logger


class ProgressBar:
    """Progress bar widget for logging.

    Parameters
    ----------
    n : int
        Total number of steps.

    logger : logging.Logger
        The `Logger` object to which to report progress.

    Examples
    --------
    with ProgressBar(1000, logger) as bar:
      for i in range(1000):
        bar.update()

    """

    def __init__(self, n, logger):
        self.n = n
        self.logger = logger

    def __enter__(self):
        self.i = 0
        self.last_tenths = 0
        self.logger.info('-' * 10)
        return self

    def __exit__(self, *args):
        print()
        self.logger.info('#' * 10)

    def update(self):
        self.i += 1
        tenths = int(self.i / self.n * 10)
        if tenths != self.last_tenths:
            self.last_tenths = tenths
            self.logger.info('#' * tenths + '-' * (10 - tenths))


class ProgressTriangle:
    """Progress triangle widget for logging.

    Parameters
    ----------
    n : int
        Total number of steps per dot.

    logger : logging.Logger
        The `Logger` object to which to report progress.

    log : bool, optional
        Use base-2 logarithmic steps.

    Examples
    --------
    with ProgressTriangle(1000, logger) as tri:
        for i in range(1000):
            tri.update()

    tri = ProgressTriangle(1, logger, log=True)
    tri.update()

    """

    def __init__(self, n, logger, log=False):
        self.n = n
        self.logger = logger
        self.log = log
        self.reset()

    def __enter__(self):
        self.reset()
        return self

    def __exit__(self, *args):
        print()

    def reset(self):
        self.i = 0

    def update(self):
        self.i += 1
        if self.log:
            logi = np.log2(self.i)
            if logi % self.n == 0:
                self.logger.info('{}'.format('.' * int(logi)))
        else:
            if self.i % self.n == 0:
                self.logger.info('{}'.format('.' * self.i))
