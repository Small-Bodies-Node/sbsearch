# Licensed with the 3-clause BSD license.  See LICENSE for details.
import sqlite3

import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord

from . import logging, util, schema
from .db import SBDB
from .config import Config


class SBSearch:
    """Search and manage survey data and small Solar System object detections.

    Parameters
    ----------
    config : sbsearch.config.Config, optional
        Use this configuration set.

    savelog : bool, optional
        Set to ``True`` to write the log to the log file.

    obs_table : string, optional
        Name of the observation table.

    obs_columns : list, optional
        Use these SQLite3 column definitions when creating the
        observation table.

    **kwargs
        If ``config`` is ``None``, pass these additional keyword
        arguments to ``Config`` initialization.

    """

    def __init__(self, config=None, savelog=False, obs_table=None,
                 obs_columns=None, **kwargs):
        self.config = Config(**kwargs) if config is None else config

        fn = self.config['log'] if savelog else '/dev/null'
        self.logger = logging.setup(filename=fn)

        self.db = sqlite3.connect(self.config['database'], 5, 0, "DEFERRED",
                                  True, SBDB)
        if obs_table is not None:
            self.db.OBS_TABLE = obs_table
        if obs_columns is not None:
            self.db.OBS_COLUMNS = obs_columns

        self.db.verify_tables(self.logger)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.logger.info('Closing database.')
        self.db.commit()
        self.db.execute('PRAGMA optimize')
        self.db.close()
        self.logger.info(Time.now().iso + 'Z')

    def add_found(self, *args, **kwargs):
        return self.db.add_found(*args, **kwargs)
    add_found.__doc__ = SBDB.add_found.__doc__

    def available_objects(self):
        """All available objects.

        Returns
        -------
        desg : ndarray
            All available object designations.

        """
        objid, desg = self.db.get_objects()
        return desg

    def clean_ephemeris(self, objects, start=None, stop=None):
        """Remove ephemeris between dates(inclusive).

        Parameters
        ----------
        objects : list
            Object IDs or designations.

        start, stop : float or `~astropy.time.Time`, optional
            Date range (inclusive), or ``None`` for unbounded.

        """
        jd_start, jd_stop = util.epochs_to_jd([start, stop])
        total = 0
        self.logger.debug('Removing ephemeris rows:')
        for objid, desg in zip([self.resolve_object(obj) for obj in objects]):
            n = self.db.clean_ephemeris(objid, jd_start, jd_stop)
            self.logger.debug('* {}: {}'.format(desg, n))
            total += n
        self.logger.info('{} ephemeris rows removed.'.format(total))

    def find_and_add(self, objects, **kwargs):
        """Find observations covering an object and add to found database.

        Parameters
        ----------
        objects : list
            Objects to find, designations and/or object IDs.

        **kwargs
           Any ``find_observations`` keyword arguments.

        """

        self.logger.info('Searching for {} objects.'.format(len(objects)))
        tri = ProgressTriangle(1, self.logger, log=True)
        for objid, desg in [self.db.resolve_object(obj) for obj in objects]:
            obsids = self.find_observations(objid, **kwargs)
            foundids = self.add_found(objid, obsids, self.config['location'])

            tri.update(len(foundids))
            self.logger.debug('* {} x{}'.format(desg, len(foundids)))

        self.logger.info('Found in {} observations.'.format(tri.i))

    def find_observations(self, obj, start=None, stop=None, vmax=22):
        """Find observations covering object.

        Parameters
        ----------
        obj : int or desg
            Object to find.

        start : float or `~astropy.time.Time`, optional
            Search after this date, inclusive.

        stop : float or `~astropy.time.Time`, optional
            Search before this date, inclusive.

        vmax : float, optional
            Require epochs brighter than this limit.

        Returns
        -------
        obsids : array of int
            Observations with this object.

        """

        objid, desg = self.db.resolve_object(obj)
        columns = ('obsid,jd_start,jd_stop,ra1,dec1,ra2,dec2,'
                   'ra3,dec3,ra4,dec4')

        segments = self.db.get_ephemeris_segments(
            objid=objid, start=start, stop=stop)

        found = []
        for ephid, segment in segments:
            obsids = self.db.get_observations_overlapping(box=segment)
            matched = self.db.get_observations_by_id(obsids, columns=columns)

            for obs in matched:
                epochs = [(obs['jd_start'] + obs['jd_stop']) / 2]
                eph, vmag = self.db.get_ephemeris_interp(objid, epochs)
                if vmag > vmax:
                    continue

                ra = np.array(obs[3::2])
                dec = np.array(obs[4::2])
                corners = SkyCoord(ra, dec, unit='rad')
                if util.interior_test(eph, corners):
                    found.append(obs['obsid'])

        return found

    def find_one_shot(self, desg, dates):
        """Find observations covering object on-the-fly.

        Parameters
        ----------
        desg : string
            Object designation.

        """
        pass

    def find_in_obs(self, obsid, exact=True):
        """Find all objects in an observation."""
        pass

    def object_coverage(self, objids, cov_type, start=None, stop=None,
                        length=60):
        """Ephemeris or found coverage for objects.

        Parameters
        ----------
        objids : list
            Object IDs to analyze.

        cov_type : string
            'eph' for ephemeris coverage, 'found' for found coverage.

        start, stop : string or astropy.time.Time, optional
            Start and stop epochs, parseable by `~astropy.time.Time`.
            Previously defined ephemerides will be removed.

        length : int, optional
            Length of the strings to create.

        Returns
        -------
        coverage : dict
            Object designations and their ephemeris coverage over the
            time period as a string:

                # = full coverage
                / = partial coverage
                - = no coverage

        """
        coverage = {}
        jd_start, jd_stop = util.epochs_to_jd(start, stop)

        rows = self.db.get_observations_by_date(
            start, stop, columns='(jd_start + jd_stop) / 2', generator=True)
        jd_obs = np.unique(list([row[0] for row in rows]))
        n_obs, bins = np.histogram(jd_obs, bins=length)

        for objid in objids:
            desg = self.db.resolve_object(objid)[1]
            if type == 'eph':
                rows = self.db.get_ephemeris(objid, jd_start, jd_stop,
                                             columns='obsjd', generator=True)
                jd = np.array(list([row[0] for row in rows]))
            elif type == 'found':
                raise NotImplemented

            ratio = np.histogram(jd, bins=bins)[0] / n_obs
            s = ''
            for r in ratio:
                if r == 0:
                    s += '-'
                elif r < 1.0:
                    s += '/'
                else:
                    s += '#'
            coverage[desg] = s

        return coverage

    def update_ephemeris(self, obj, start, stop, step=None, cache=False):
        """Update object ephemeris table.

        Parameters
        ----------
        obj : string or integer
            Designation resolvable by the Minor Planet Center, or
            object table ID.

        start, stop : string or astropy.time.Time
            Start and stop epochs, parseable by `~astropy.time.Time`.
            Previously defined ephemerides will be removed.

        step : string or astropy.unit.Quantity, optional
            Integer step size in hours or days.  If `None`, then an
            adaptable step size will be used, based on distance to the
            observer.

        cache : bool, optional
            ``True`` to use ``astroquery`` cache, primarily for
            testing.

        """

        objid = self.db.resolve_object(obj)[0]
        jd_start, jd_stop = util.epochs_to_jd([start, stop])

        n = self.db.clean_ephemeris(objid, jd_start, jd_stop)
        self.logger.info('{} rows deleted from eph table'.format(n))

        n = self.db.add_ephemeris(
            objid, self.config['location'], jd_start, jd_stop, step=step)

        self.logger.info('{} rows added to eph table'.format(n))
