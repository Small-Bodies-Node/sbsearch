# Licensed under a 3-clause BSD style license - see LICENSE.rst
import abc
import sys
import logging
from typing import Callable, Optional
import numpy as np
from astropy.time import Time


def setup_logger(filename: str = 'sbsearch.log', name: str = 'SBSearch',
                 level: Optional[int] = None) -> logging.Logger:
    logger: logging.Logger = logging.Logger(name)
    if level:
        logger.setLevel(level)
    else:
        logger.setLevel(logging.DEBUG)

    # delete any previous handlers
    logger.handlers = []

    formatter: logging.Formatter = logging.Formatter(
        '%(levelname)s %(asctime)s (' + name + '): %(message)s')

    console: logging.StreamHandler = logging.StreamHandler(sys.stdout)
    if level is None:
        console.setLevel(logging.DEBUG)
    console.setFormatter(formatter)
    logger.addHandler(console)

    logfile: logging.FileHandler = logging.FileHandler(filename)
    if level is None:
        logfile.setLevel(logging.INFO)
    logfile.setFormatter(formatter)
    logger.addHandler(logfile)

    return logger


class ProgressWidget(abc.ABC):
    """Log progress as percentage.


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

    def __init__(self, n: int, logger: logging.Logger) -> None:
        self.n: int = n
        self.logger: logging.Logger = logger
        self.i: int = 0
        self.t0: Time = Time.now()

    def __enter__(self):
        self.reset()
        return self

    def __exit__(self, *args):
        self.done()

    @ property
    def dt(self):
        return (Time.now() - self.t0).sec

    def reset(self):
        """Reset the counter."""
        self.i = 0
        self.t0: Time = Time.now()

    def log(self, status: float):
        """Write the current status to the log."""
        self.logger.info('%.2f%%', status * 100)

    def update(self, increment: int = 1):
        """Increment counter and report status."""
        self.i += increment
        self.log(self.i / self.n)

    def done(self):
        self.logger.info('{:.0f} seconds elapsed.'.format(self.dt))


class ProgressBar(ProgressWidget):
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

    def __init__(self, n: int, logger: logging.Logger) -> None:
        super().__init__(n, logger)
        self.last_tenths: int = 0

    def reset(self):
        self.last_tenths = 0
        self.log(0)
        super().reset()

    def log(self, tenths: int) -> None:
        self.logger.info('%s', '#' * tenths + '-' * (10 - tenths))

    def update(self, increment: int = 1):
        self.i += increment
        tenths: int = int(self.i / self.n * 10)
        if tenths != self.last_tenths:
            self.last_tenths = tenths
            self.log(tenths)


class ProgressTriangle(ProgressWidget):
    """Progress triangle widget for logging.


    Parameters
    ----------
    n : int
        Total number of steps per dot.

    logger : logging.Logger
        The `Logger` object to which to report progress.

    base : int, optional
        Use logarithmic steps with this base: 2 or 10.


    Examples
    --------
    with ProgressTriangle(1000, logger) as tri:
        for i in range(1000):
            tri.update()

    """

    def __init__(self, n: int, logger: logging.Logger, base: int = 2) -> None:
        super().__init__(n, logger)
        self.base = base
        self._log: Callable
        if base == 2:
            self._log = np.log2
            self.logger.info('Base-2 dots')
        elif base == 10:
            self._log = np.log10
            self.logger.info('Base-10 dots')
        else:
            raise ValueError('base must be 2 or 10')

    def log(self, status: int) -> None:
        self.logger.info('%s', '.' * status)

    def update(self, increment: int = 1):
        last = self.i
        self.i += increment
        if last == 0:
            return

        logi = self._log(self.i)
        if (self._log(last) % self.n) >= (logi % self.n):
            self.log(int(logi // self.n))
