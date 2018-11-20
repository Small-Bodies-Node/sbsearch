# Licensed under a 3-clause BSD style license - see LICENSE.rst

import logging as pylogging
import pytest
from .. import logging


def test_setup():
    logger = logging.setup()
    assert logger.level == pylogging.DEBUG

    logger = logging.setup(level=pylogging.ERROR)
    assert logger.level == pylogging.ERROR


def test_ProgressBar_stdout(capsys):
    N = 137
    with logging.ProgressBar(N) as progress:
        for i in range(N):
            progress.update()
    captured = capsys.readouterr()
    assert captured.out.strip() == '''----------
#---------
##--------
###-------
####------
#####-----
######----
#######---
########--
#########-
##########
'''.strip()


def test_ProgressBar_log(caplog):
    N = 137
    caplog.set_level(pylogging.INFO)
    logger = pylogging.getLogger('progress')
    with logging.ProgressBar(N, logger=logger) as progress:
        for i in range(N):
            progress.update()

    assert caplog.record_tuples == [
        ('progress', 20, '----------'), ('progress', 20, '#---------'),
        ('progress', 20, '##--------'), ('progress', 20, '###-------'),
        ('progress', 20, '####------'), ('progress', 20, '#####-----'),
        ('progress', 20, '######----'), ('progress', 20, '#######---'),
        ('progress', 20, '########--'), ('progress', 20, '#########-'),
        ('progress', 20, '##########')]


def test_ProgressTriangle_stdout_linear(capsys):
    N = 100
    with logging.ProgressTriangle(20) as progress:
        for i in range(N):
            progress.update()
    captured = capsys.readouterr()

    # avoid comparing timing column
    lines = captured.out.strip().splitlines()[:-1]
    dots = [line.split()[1] for line in lines]
    assert dots == '''.
..
...
....
.....'''.splitlines()


def test_ProgressTriangle_logger_log2(caplog):
    N = 100
    caplog.set_level(pylogging.INFO)
    logger = pylogging.getLogger('progress')
    with logging.ProgressTriangle(1, logger=logger, base=2) as progress:
        for i in range(N):
            progress.update()

    # avoid comparing timing column
    expected = '''.
..
...
....
.....
......'''.splitlines()
    for record, test in zip(caplog.record_tuples[1:-1], expected):
        assert record[2].split()[-1] == test


def test_ProgressTriangle_logger_log10(caplog):
    N = 1000
    caplog.set_level(pylogging.INFO)
    logger = pylogging.getLogger('progress')
    progress = logging.ProgressTriangle(1, logger=logger, base=10)
    for i in range(N):
        progress.update()
    progress.done()
    # avoid comparing timing column
    expected = '''.
..
...'''.splitlines()
    for record, test in zip(caplog.record_tuples[1:-1], expected):
        assert record[2].split()[-1] == test


def test_ProgressTriangle_error():
    with pytest.raises(ValueError):
        progress = logging.ProgressTriangle(1, base=3)
