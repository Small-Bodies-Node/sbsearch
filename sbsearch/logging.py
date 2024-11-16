# Licensed under a 3-clause BSD style license - see LICENSE.rst
import abc
import sys
import time
import logging
from enum import Enum
from typing import Callable, Optional, Union
import numpy as np
from astropy.time import Time


class ElapsedFormatter(logging.Formatter):
    """Based on https://stackoverflow.com/questions/25194864/python-logging-time-since-start-of-program"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.t0 = time.time()
        self.t_last = self.t0

    def format(self, record):
        record.dt0 = record.created - self.t0
        record.dt = record.created - self.t_last
        self.t_last = record.created
        return super().format(record)


def setup_logger(
    filename: str = "sbsearch.log", name: str = "SBSearch", level: Optional[int] = None
) -> logging.Logger:
    logger: logging.Logger = logging.getLogger(name)

    if level:
        logger.setLevel(level)
    else:
        logger.setLevel(logging.DEBUG)

    # reset handlers, in case already defined
    close_logger(name)

    formatter = ElapsedFormatter(
        "%(asctime)10s (%(dt).2f/%(dt0).2f) %(levelname)s: "
        "[%(funcName)s] %(message)s",
        "%Y-%m-%d %H:%M:%S",
    )

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


def close_logger(name: str = "SBSearch") -> None:
    """Close SBSearch loggers."""
    logger: logging.Logger = logging.getLogger(name)
    for handler in logger.handlers:
        handler.close()
        logger.removeHandler(handler)


class ProgressWidget(abc.ABC):
    def __init__(self, n: int, logger: logging.Logger) -> None:
        self.n: int = n
        self.logger: logging.Logger = logger
        self.i: int = 0
        self.t0: float = time.monotonic()

    def __enter__(self):
        self.reset()
        return self

    def __exit__(self, *args):
        self.done()

    @property
    def dt(self):
        return int(time.monotonic() - self.t0)

    def reset(self):
        """Reset the counter."""
        self.i = 0
        self.t0 = time.monotonic()

    @abc.abstractmethod
    def log(self):
        pass

    @abc.abstractmethod
    def update(self):
        pass

    def done(self):
        self.logger.info("%d in %.0f seconds", self.i, self.dt)


class ProgressPercent(ProgressWidget):
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

    def log(self, status: float):
        """Write the current status to the log."""
        self.logger.info("%.2f%%", status * 100)

    def update(self, increment: int = 1):
        """Increment counter and report status."""
        self.i += increment
        self.log(self.i / self.n)


class ProgressScale(Enum):
    LOG = "log"
    LINEAR = "linear"


class ProgressBar(ProgressWidget):
    """Progress bar widget for logging.


    Parameters
    ----------
    n : int
        Total number of steps.

    logger : logging.Logger
        The `Logger` object to which to report progress.

    scale : string, optional
        'log' or 'linear' scale.

    length : int, optional
        Total length of the progress bar.


    Examples
    --------
    with ProgressBar(1000, logger) as bar:
        for i in range(1000):
            bar.update()

    """

    def __init__(
        self, n: int, logger: logging.Logger, scale: str = "linear", length: int = 10
    ) -> None:
        super().__init__(n, logger)
        self.last_step: int = 0
        self.scale: ProgressScale = ProgressScale(scale)
        self.length: int = length

    def reset(self):
        super().reset()
        self.last_step = 0
        self.logger.info("%d steps, %s scale", self.length, self.scale.value)
        self.logger.info("%s", "-" * self.length)

    def log(self, step: int) -> None:
        self.logger.info(
            "%s t-%d s",
            "#" * step + "-" * (self.length - step),
            int((time.monotonic() - self.t0) / self.i * (self.n - self.i)),
        )

    def update(self, increment: int = 1):
        self.i += increment

        step: int
        if self.scale == ProgressScale.LINEAR:
            step = int(self.i / self.n * self.length)
        else:
            try:
                step = int(np.log10(self.i) / np.log10(self.n) * self.length)
            except ValueError:
                if self.i == self.n:
                    step = self.length
                else:
                    step = 0

        if step != self.last_step:
            self.last_step = step
            self.log(step)


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
            self.logger.info("Base-2 dots")
        elif base == 10:
            self._log = np.log10
            self.logger.info("Base-10 dots")
        else:
            raise ValueError("base must be 2 or 10")

    def log(self, status: int) -> None:
        self.logger.info("%s %d", "." * status, self.i)

    def update(self, increment: int = 1):
        last = self.i
        self.i += increment
        if last == 0:
            return

        logi = self._log(self.i)
        if (self._log(last) % self.n) >= (logi % self.n):
            self.log(int(logi // self.n))
