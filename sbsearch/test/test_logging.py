# Licensed under a 3-clause BSD style license - see LICENSE.rst

import logging
import pytest
from ..logging import(setup_logger, close_logger, ProgressBar,
                      ProgressTriangle, ProgressWidget)


def test_setup():
    logger = setup_logger()
    assert logger.level == logging.DEBUG
    close_logger()
    
    logger = setup_logger(level=logging.ERROR)
    assert logger.level == logging.ERROR
    close_logger()


def test_ProgressWidget_log(caplog):
    N = 6
    caplog.set_level(logging.INFO)
    logger = logging.getLogger('progress')
    with ProgressWidget(N, logger) as progress:
        for i in range(N):
            progress.update()

    assert caplog.record_tuples == [
        ('progress', 20, '16.67%'), ('progress', 20, '33.33%'),
        ('progress', 20, '50.00%'), ('progress', 20, '66.67%'),
        ('progress', 20, '83.33%'), ('progress', 20, '100.00%'),
        ('progress', 20, '0 seconds elapsed.')
    ]


def test_ProgressBar_log(caplog):
    N = 137
    caplog.set_level(logging.INFO)
    logger = logging.getLogger('progress')
    with ProgressBar(N, logger=logger) as progress:
        for i in range(N):
            progress.update()

    assert caplog.record_tuples == [
        ('progress', 20, '10 steps, linear scale'),
        ('progress', 20, '----------'),
        ('progress', 20, '#--------- t-0 s'),
        ('progress', 20, '##-------- t-0 s'),
        ('progress', 20, '###------- t-0 s'),
        ('progress', 20, '####------ t-0 s'),
        ('progress', 20, '#####----- t-0 s'),
        ('progress', 20, '######---- t-0 s'),
        ('progress', 20, '#######--- t-0 s'),
        ('progress', 20, '########-- t-0 s'),
        ('progress', 20, '#########- t-0 s'),
        ('progress', 20, '########## t-0 s'),
        ('progress', 20, '0 seconds elapsed.')]


def test_ProgressTriangle_logger_log2(caplog):
    N = 100
    caplog.set_level(logging.INFO)
    logger = logging.getLogger('progress')
    with ProgressTriangle(1, logger, base=2) as progress:
        for i in range(N):
            progress.update()

    expected = '''. 2
.. 4
... 8
.... 16
..... 32
...... 64'''.splitlines()
    for record, test in zip(caplog.record_tuples[1:-1], expected):
        assert record[2].strip() == test


def test_ProgressTriangle_logger_log10(caplog):
    N = 1000
    caplog.set_level(logging.INFO)
    logger = logging.getLogger('progress')
    progress = ProgressTriangle(1, logger, base=10)
    for i in range(N):
        progress.update()
    progress.done()
    # avoid comparing timing column
    expected = '''. 10
.. 100
... 1000'''.splitlines()
    for record, test in zip(caplog.record_tuples[1:-1], expected):
        assert record[2].strip() == test


def test_ProgressTriangle_error():
    logger = logging.getLogger('progress')
    with pytest.raises(ValueError):
        ProgressTriangle(1, logger, base=3)
