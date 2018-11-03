# Licensed with the 3-clause BSD license.  See LICENSE for details.
import re
import itertools

import numpy as np
import astropy.units as u
from astropy.coorindates import SkyCoord

from . import logging, util, schema
from .config import Config
from .execptions import CorruptDB


class SBSearch:
    """Search and manage survey data and small Solar System object detections.

    Parameters
    ----------
    config : sbsearch.config.Config, optional
        Use this configuration set.

    savelog : bool, optional
        Set to ``True`` to write the log to the log file.

    """

    def __init__(self, config=None, savelog=False):
        self.config = Config() if config is None else config
        fn = self.config['logfile'] if log else '/dev/null'
        self.logger = logging.setup(filename=fn)
        self.db = sqilte3.connect(config['filename'], factory=SBDB)
        self.db.verify_tables(self.config, self.logger)

    def update_obs(self, observations):
        """Update observation table.

        Parameters
        ----------
        observations : list of array-like
            List of rows to insert (or replace) into the observation
            database.  All database columns must have a corresponding
            value and the values must match the table column order.

        """
        count = 0
        for observation in observations:
            c = self.db.executemany('''
            INSERT OR REPLACE INTO ? VALUES ({})
            '''.format(','.join('?' * len(observation))),
                [self.config['obs_table']] + list(observation))
            count += c.rowcount

        self.logger.info(
            'Updated observation table with {} rows'.format(count))

    def update_eph(self, obj, start, stop, step=None, clean=False):
        """Update object ephemeris table.

        Parameters
        ----------
        obj : string or integer
            Designation resolvable by the Minor Planet Center, or
            object table ID.

        start, stop : string or astropy.time.Time
            Start and stop epochs, parseable by `~astropy.time.Time`.

        step : string or astropy.unit.Quantity, optional
            Integer step size in hours or days.  If `None`, then an
            adaptable step size will be used, based on distance to the
            Earth.

        clean : bool, optional
            If ``True``, remove ephemerides, between and including
            start and stop, from the database.  With ``start is
            None``, remove all points before ``stop``.  With ``stop is
            None``, remove all points after ``start``.  If both are
            ``None``, all ephemeris points are removed.

        """

        objid = self.db.resolve_object(obj)[0]
        jd_start = None if start is None else Time(start, scale='utc').jd
        jd_stop = None if stop is None else Time(stop, scale='utc').jd
        if clean:
            n = self.db.clean_ephemeris(objid, jd_start, jd_stop)
            self.logger.info('{} rows deleted from eph table'.format(n))
        else:
            n = self.db.add_ephemeris(
                self.config['location'], objid, jd_start, jd_stop, step=step)
            self.logger.info('{} rows added to eph table'.format(n))

    def clean_found(self, obj, start=None, stop=None):
        """Remove found objects from the database.

        Parameters
        ----------
        obj : str or int
            Object to remove as a designation or obsid.

        start, stop : string or astropy.time.Time, optional
            Remove epochs between these times, unbounded if ``None``.

        """

        objid = self.resolve_object(obj)[0]
        jd_start = None if start is None else Time(start, scale='utc').jd
        jd_stop = None if stop is None else Time(stop, scale='utc').jd

        constraints = [('objid=?', objid)]
        constraints.extend(util.date_constraints(jd_start, jd_stop))
        cmd, parameters = util.assemble_sql(
            'DELETE FROM ' + self.config.found_table, [], constaints)

        c = self.db.execute(cmd, parameters)
        self.logger.info('{} rows deleted from {}'.format(
            c.rowcount, self.config.found_table))

    def find_by_date(self, desg, dates, exact=True):
        """Find observations covering objects."""
        pass

    def find_in_obs(self, obsid, exact=True):
        """Find all objects in an observation."""
        pass
