# Licensed with the 3-clause BSD license.  See LICENSE for details.
import sqlite3
from itertools import repeat

import numpy as np
from astropy.time import Time
from astropy.table import Table, Column, vstack

from . import logging, util, schema
from .util import RADec
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

    schema : str, optional
        Additional database schema.

    **kwargs
        If ``config`` is ``None``, pass these additional keyword
        arguments to ``Config`` initialization.

    """

    def __init__(self, config=None, savelog=False, **kwargs):
        self.config = Config(**kwargs) if config is None else config
        self.config.update(**kwargs)

        fn = self.config['log'] if savelog else '/dev/null'
        self.logger = logging.setup(filename=fn)

        self.db = sqlite3.connect(self.config['database'], 5, 0, "DEFERRED",
                                  True, SBDB)
        self.verify_database()

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
        """Remove ephemeris between dates (inclusive).

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
        for objid, desg in self.db.resolve_objects(objects):
            if objid is None:
                continue
            n = self.db.clean_ephemeris(objid, jd_start, jd_stop)
            self.logger.debug('* {}: {}'.format(desg, n))
            total += n
        self.logger.info('{} ephemeris rows removed.'.format(total))

    def clean_found(self, objects, start=None, stop=None):
        """Remove found entries between dates (inclusive).

        Parameters
        ----------
        objects : list
            Object IDs or designations, or ``None`` to remove all
            found entries.

        start, stop : float or `~astropy.time.Time`, optional
            Date range (inclusive), or ``None`` for unbounded.

        """
        jd_start, jd_stop = util.epochs_to_jd([start, stop])
        total = 0
        self.logger.debug('Removing found rows:')
        if objects is None:
            total = self.db.clean_found(objid, jd_start, jd_stop)
        else:
            for objid, desg in self.db.resolve_objects(objects):
                if objid is None:
                    self.logger.warning(
                        '{} not in object table'.format(desg))
                    continue

                n = self.db.clean_found(objid, jd_start, jd_stop)
                self.logger.debug('* {}: {}'.format(desg, n))
                total += n

        self.logger.info('{} found rows removed.'.format(total))

    def find_objects(self, objects, start=None, stop=None, vmax=25,
                     save=False, update=False):
        """Find observations covering an object.

        Tests for observations in requested range before searching.

        Parameters
        ----------
        objects : list or ``None``
            Objects to find, designations and/or object IDs.  ``None``
            to search for all objects.

        start : float or `~astropy.time.Time`, optional
            Search after this epoch, inclusive.

        stop : float or `~astropy.time.Time`, optional
            Search before this epoch, inclusive.

        vmax : float, optional
            Require epochs brighter than this limit.

        save : bool, optional
            If ``True``, store observations to local database.

        update : bool, optional
            If ``True``, update metadata for already found objects.

        Returns
        -------
        tab : `~astropy.table.Table`
            Summary of found objects.

        """

        if not self.db.test_observation_coverage(start, stop):
            self.logger.info(
                'No observations in database over requested range.')
            return self.observation_summary([])

        if start is None and stop is None:
            s = 'in all observations'
        elif start is None:
            d = util.epochs_to_time([stop])[0].iso
            s = 'in observations ending {}'.format(d)
        elif stop is None:
            d = util.epochs_to_time([start])[0].iso
            s = 'in observations starting {}'.format(d)
        elif start == stop:
            d = util.epochs_to_time([start])[0].iso[:10]
            s = 'in observations on {}'.format(d)
        else:
            df, dt = util.epochs_to_time([start, stop]).iso
            s = 'in observations from {} to {}'.format(df, dt)

        s += ', V<={:.1f}'.format(vmax)

        self.logger.info('Searching for {} object{} {}.'.format(
            len(objects), 's' if len(objects) == 1 else '', s))

        n = 0
        summary = []
        progress = logging.ProgressTriangle(1, self.logger, base=2)
        for objid, desg in self.db.resolve_objects(objects):
            obs = self.find_object(objid, start=start, stop=stop,
                                   vmax=vmax)
            if len(obs) == 0:
                continue

            n += len(obs)

            if save:
                obsids = [o['obsid'] for o in obs]
                foundids = self.add_found(
                    objid, obsids, self.config['location'], update=update)

            tab = self.observation_summary(obsids, add_found=save)
            tab.add_column(Column(list(repeat(desg, len(tab))), name='desg'),
                           index=0)
            summary.append(tab)

            progress.update(len(obs))
            self.logger.debug('* {} x{} ({} saved)'.format(
                desg, len(obs), len(foundids)))

        progress.done()
        self.logger.info(
            'Found in {} observations ({} searched in detail).'.format(
                progress.i, n))

        if len(summary) == 0:
            return self.observation_summary([])
        else:
            return vstack(summary)

    def find_object(self, obj, start=None, stop=None, vmax=25,
                    source=None):
        """Find observations covering object.

        Does not test for observation coverage before search.

        Parameters
        ----------
        obj : int or desg
            Object to find.

        start : float or `~astropy.time.Time`, optional
            Search after this epoch, inclusive.

        stop : float or `~astropy.time.Time`, optional
            Search before this epoch, inclusive.

        vmax : float, optional
            Require epochs brighter than this limit.

        source : string, optional
            Ephemeris source: ``None`` for internal database, 'mpc' or
            'jpl' for online ephemeris generation.

        Returns
        -------
        obs : tuple
            Observations with this object.

        """

        objid, desg = self.db.resolve_object(obj)
        start, stop = util.epochs_to_jd((start, stop))

        segments = self.db.get_ephemeris_segments(
            objid=objid, start=start, stop=stop)

        found = []
        for ephid, segment in segments:
            obsids = self.db.get_observations_overlapping(box=segment)
            matched = self.db.get_observations_by_id(obsids)

            for obs in matched:
                epoch = [(obs['jd_start'] + obs['jd_stop']) / 2]
                eph, vmag = self.db.get_ephemeris_interp(objid, epoch)
                if vmag > vmax:
                    continue

                ra = np.array([obs['ra' + i] for i in '1234'])
                dec = np.array([obs['dec' + i] for i in '1234'])
                point = RADec(eph.ra, eph.dec, unit='rad')
                corners = RADec(ra, dec, unit='rad')
                if util.interior_test(point, corners):
                    found.append(obs)

        return list(set(found))

    def find_by_ephemeris(self, eph):
        """Find object based on ephemeris.

        Use for objects that are not in the database.

        Parameters
        ----------
        eph : `~sbpy.data.Ephem`
            Ephemeris of object to find, requires RA, Dec, and date columns.

        Returns
        -------
        n : int
            Number of observations checked in detail.

        obs : tuple
            Observations with this object.

        tab : `~astropy.table.Table`
            Summary of found observations.

        """

        found = []
        n = 0
        jd = np.array(util.epochs_to_jd(eph['Date'].value))
        mjd = jd - 2400000.5
        coords = RADec(eph['RA'], eph['Dec'])
        x, y, z = coords.xyz
        for i in range(len(eph) - 1):
            segment = {
                'mjd0': mjd[i],
                'mjd1': mjd[i + 1],
                'x0': min(x[i:i+2]),
                'x1': max(x[i:i+2]),
                'y0': min(y[i:i+2]),
                'y1': max(y[i:i+2]),
                'z0': min(z[i:i+2]),
                'z1': max(z[i:i+2])
            }

            obsids = self.db.get_observations_overlapping(box=segment)
            matched = self.db.get_observations_by_id(obsids)

            for obs in matched:
                jd0 = (obs['jd_start'] + obs['jd_stop']) / 2

                point = util.spherical_interpolation(
                    coords[i], coords[i + 1], jd[i], jd[i + 1], jd0)

                ra = np.array([obs['ra' + i] for i in '1234'])
                dec = np.array([obs['dec' + i] for i in '1234'])
                corners = RADec(ra, dec, unit='rad')
                if util.interior_test(point, corners):
                    found.append(obs)
                n += 1

        tab = self.observation_summary(found)

        self.logger.info('{} ephemeris positions retrieved.'.format(len(eph)))
        self.logger.info('{} observations searched in detail.'.format(n))
        self.logger.info('{} observations found.'.format(len(found)))
        return n, found, tab

    def find_in_observation(self, obsid, exact=True):
        """Find all objects in a single observation."""
        pass

    def object_coverage(self, objids, cov_type, start=None, stop=None,
                        length=60):
        """Ephemeris or found coverage for objects.

        Parameters
        ----------
        objids: list
            Object IDs to analyze.

        cov_type: string
            'eph' for ephemeris coverage, 'found' for found coverage.

        start, stop: string or astropy.time.Time, optional
            Start and stop epochs, parseable by `~astropy.time.Time`.
            Previously defined ephemerides will be removed.

        length: int, optional
            Length of the strings to create.

        Returns
        -------
        coverage: dict
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

    def observation_summary(self, obsids, add_found=False):
        """Summarize observations.

        Parameters
        ----------
        obsids : tuple or list of int
            Observation IDs to summarize.

        add_found : bool, optional
            Add metadata from found table.

        Returns
        -------
        summary : Table
            Summary table.

        """

        names = ('obsid', 'date', 'ra', 'dec')
        if add_found:
            columns = ('obsid,(jd_start + jd_stop) / 2 AS jd,'
                       'found.ra,found.dec,rh,delta,vmag')
            names += ('rh', 'delta', 'vmag')
            inner_join = ['found USING obsid']
        else:
            columns = 'obsid,(jd_start + jd_stop) / 2 AS jd,ra,dec'
            inner_join = None

        rows = []
        obs = self.db.get_observations_by_id(
            obsids, inner_join=inner_join, generator=True)
        for row in obs:
            rows.append([row['obsid'], Time(row['jd'], format='jd').iso[:-4]]
                        + list([r for r in row[2:]]))

        return Table(rows=rows, names=names)

    def update_ephemeris(self, objects, start, stop, step=None, cache=False):
        """Update object ephemeris table.

        Parameters
        ----------
        objects: list of string or integer
            Designations resolvable by the Minor Planet Center, or
            object table IDs.

        start, stop: string or astropy.time.Time
            Start and stop epochs, parseable by `~astropy.time.Time`.
            Previously defined ephemerides will be removed.

        step: string or astropy.unit.Quantity, optional
            Integer step size in hours or days.  If `None`, then an
            adaptable step size will be used, based on distance to the
            observer.

        cache: bool, optional
            ``True`` to use ``astroquery`` cache, primarily for
            testing.

        """

        jd_start, jd_stop = util.epochs_to_jd([start, stop])

        cleaned = 0
        added = 0
        for objid, desg in self.db.resolve_objects(objects):
            self.logger.debug('* {}'.format(desg))
            if objid is None:
                objid = self.db.add_object(desg)

            n = self.db.clean_ephemeris(objid, jd_start, jd_stop)
            cleaned += n

            n = self.db.add_ephemeris(
                objid, self.config['location'], jd_start, jd_stop, step=step)
            added += n

        self.logger.info('{} rows deleted from eph table'.format(cleaned))
        self.logger.info('{} rows added to eph table'.format(added))

    def verify_database(self, names=[], script=''):
        """Verify database tables, triggers, etc."""
        self.db.verify_database(self.logger, names=names, script=script)
