# Licensed with the 3-clause BSD license.  See LICENSE for details.
import sqlite3
from itertools import repeat
from logging import ERROR
import requests

import numpy as np
from astropy.io import ascii
from astropy.time import Time
from astropy.table import Table, Column, vstack
import astropy.units as u
from sbpy.data import Orbit, Ephem

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

    save_log : bool, optional
        Set to ``True`` to write the log to the log file.

    disable_log : bool, optional
        Set to ``True`` to disable normal logging; also sets
        ``save_log=False``.

    **kwargs
        If ``config`` is ``None``, pass these additional keyword
        arguments to ``Config`` initialization.

    """

    def __init__(self, config=None, save_log=False, disable_log=False,
                 **kwargs):
        self.config = Config(**kwargs) if config is None else config
        self.config.update(**kwargs)

        if disable_log:
            save_log = False
            level = ERROR
        else:
            level = None

        fn = self.config['log'] if save_log else '/dev/null'
        self.logger = logging.setup(filename=fn, level=level)

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
        """List available objects.

        Returns
        -------
        tab : `~astropy.table.Table`

        """
        tab = Table(self.db.get_objects(),
                    names=('object ID', 'designation'))
        return tab

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
                self.logger.warning('{} not in object table'.format(desg))
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
            total = self.db.clean_found(None, jd_start, jd_stop)
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

        found : list
            Observation IDs with this object.

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

                ra, dec = util.fov2points(obs['fov'])
                corners = RADec(ra[1:], dec[1:], unit='rad')
                if util.interior_test(point, corners):
                    found.append(obs['obsid'])
                n += 1

        tab = self.observation_summary(found)

        self.logger.info('{} observations searched in detail.'.format(n))
        self.logger.info('{} observations found.'.format(len(found)))
        return n, found, tab

#    def find_in_observation(self, obsid, exact=True):
#        """Find all objects in a single observation."""
#        pass

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

                point = RADec(eph.ra, eph.dec, unit='rad')
                ra, dec = util.fov2points(obs['fov'])
                corners = RADec(ra[1:], dec[1:], unit='rad')
                if util.interior_test(point, corners):
                    found.append(obs)

        return list(set(found))

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

        if not self.db.get_observations_overlapping(start=start, stop=stop):
            self.logger.info(
                'No observations in database over requested range.')
            return self.observation_summary([])

        if start is None and stop is None:
            s = 'in all observations'
        elif start is None:
            d = util.epochs_to_time([stop])[0].iso[:16]
            s = 'in observations ending {} UT'.format(d)
        elif stop is None:
            d = util.epochs_to_time([start])[0].iso[:16]
            s = 'in observations starting {} UT'.format(d)
        elif start == stop:
            d = util.epochs_to_time([start])[0].iso[:10]
            s = 'in observations on {}'.format(d)
        else:
            df, dt = util.epochs_to_time([start, stop]).iso
            df = df[:16]
            dt = dt[:16]
            s = 'in observations from {} UT to {} UT'.format(df, dt)

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

            obsids = [o['obsid'] for o in obs]
            n += len(obs)

            if save:
                foundids = self.add_found(
                    objid, obsids, self.config['location'], update=update)
            else:
                foundids = []

            tab = self.observation_summary(obsids, add_found=save)
            tab.add_column(Column(list(repeat(desg, len(tab))), name='desg'),
                           index=0)
            summary.append(tab)

            progress.update(len(obs))
            self.logger.debug('* {} x{}, {} saved'.format(
                desg, len(obs), len(foundids)))

        progress.done()
        self.logger.info(
            'Found in {} observations ({} searched in detail).'.format(
                progress.i, n))

        if len(summary) == 0:
            return self.observation_summary([])
        else:
            return vstack(summary)

    def object_coverage(self, cov_type, objects=None, start=None, stop=None,
                        length=60):
        """Ephemeris or found coverage for objects.


        Parameters
        ----------
        cov_type: string
            'eph' for ephemeris coverage, 'found' for found coverage.

        objects: list, optional
            Object IDs or desginations to analyze.

        start, stop: string or astropy.time.Time, optional
            Start and stop epochs, parseable by `~astropy.time.Time`.

        length: int, optional
            Length of the strings to create.


        Returns
        -------
        tab : `~astropy.table.Table`

        """

        if objects is None:
            objects = zip(*self.db.get_objects())
        else:
            objects = self.db.resolve_objects(objects)

        jd_start, jd_stop = util.epochs_to_jd((start, stop))

        if jd_start is None:
            r = self.db.execute('''
            SELECT mjd0 FROM obs_tree ORDER BY mjd0 LIMIT 1
            ''').fetchone()
            jd_start = r[0] + 2400000.5

        if jd_stop is None:
            r = self.db.execute('''
            SELECT mdj1 FROM obs_tree ORDER BY mjd1 DESC LIMIT 1
            ''').fetchone()
            jd_stop = r[0] + 2400000.5

        bins = np.linspace(jd_start, jd_stop, length + 1)

        coverage = []
        for objid, desg in objects:
            if objid is None:
                coverage.append(('', desg, '-' * length))
                continue

            if cov_type == 'eph':
                rows = self.db.get_ephemeris(objid, jd_start, jd_stop,
                                             columns='jd')
            elif cov_type == 'found':
                rows = self.db.get_found(obj=objid, start=jd_start,
                                         stop=jd_stop, columns='obsjd')
            else:
                raise ValueError(
                    'coverage type must be eph or found: {}'.format(
                        cov_type))

            jd = np.array(list([row[0] for row in rows]))
            count = np.histogram(jd, bins=bins)[0]
            s = ''
            for c in count:
                if c > 0:
                    s += '+'
                else:
                    s += '-'

            coverage.append((objid, desg, s))

        tab = Table(rows=coverage, masked=True,
                    names=('object ID', 'desgination', 'coverage'))
        tab.meta['start date'] = Time(jd_start, format='jd').iso[:19]
        tab.meta['stop date'] = Time(jd_stop, format='jd').iso[:19]

        dt = (jd_stop - jd_start) / length
        x = []
        for scale, label in ((1, 'd'), (24, 'h'), (60, 'm'), (60, 's')):
            n, dt = np.divmod(dt * scale, 1)
            x.append('{}{}'.format(int(n), label))
        tab.meta['time per pip'] = ' '.join(x)
        return tab

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
            columns = ('obsid,(jd_start + jd_stop) / 2 AS jd,ra,dec,'
                       'rh,delta,vmag')
            names += ('rh', 'delta', 'vmag')
            inner_join = ['found USING obsid']
        else:
            columns = 'obsid,(jd_start + jd_stop) / 2 AS jd,fov'
            inner_join = None

        rows = []
        observations = self.db.get_observations_by_id(
            obsids, columns=columns, inner_join=inner_join, generator=True)
        for obs in observations:
            if add_found:
                ra = obs['ra']
                dec = obs['dec']
            else:
                p = util.fov2points(obs['fov'])
                ra, dec = p[0][0], p[1][0]

            rows.append([obs['obsid'],
                         Time(obs['jd'], format='jd').iso[:-4],
                         ra, dec]
                        + list([r for r in obs[4:]]))

        if len(rows) == 0:
            return Table(dtype=[float] * len(names), names=names)
        else:
            return Table(rows=rows, names=names)

    def pccp_check(self, start=None, stop=None):
        """Search for today's objects on the MPC's PCCP.

        Possible Comet Confirmation Page:
        https://minorplanetcenter.net/iau/NEO/pccp_tabular.html

        Parameters
        ----------
        start : float or `~astropy.time.Time`, optional
            Search after this epoch, inclusive.

        stop : float or `~astropy.time.Time`, optional
            Search before this epoch, inclusive.

        Returns
        -------
        tab : `~astropy.table.Table`

        """

        r = requests.get('https://minorplanetcenter.net/iau/NEO/pccp.txt')
        pccp = r.content.decode()

        if len(pccp) == 0:
            return self.observation_summary([])

        jd_start, jd_stop = util.epochs_to_jd((start, stop))

        if jd_start is None:
            r = self.db.execute('''
            SELECT mjd0 FROM obs_tree ORDER BY mjd0 LIMIT 1
            ''').fetchone()
            jd_start = r[0] + 2400000.5

        if jd_stop is None:
            r = self.db.execute('''
            SELECT mdj1 FROM obs_tree ORDER BY mjd1 DESC LIMIT 1
            ''').fetchone()
            jd_stop = r[0] + 2400000.5

        dt = jd_stop - jd_start
        steps = int(dt)
        epochs = start.jd + np.arange(steps)
        if epochs[-1] != stop.jd:
            epochs = np.r_[epochs, stop.jd]

        desgs = [line.split()[0] for line in pccp.splitlines()]
        summary = []
        for desg in desgs:
            r = requests.get(
                'https://cgi.minorplanetcenter.net/cgi-bin/showobsorbs.cgi?'
                'Obj={}&orb=y'.format(desg))
            text = r.content.decode()
            # remove html
            text = '\n'.join(text.splitlines()[1:-1])

            try:
                tab = ascii.read(text)
            except ascii.InconsistentTableError:
                self.logger.error('Cannot read orbital parameters for {}'.format(
                    desg))
                continue

            if tab['a'][0] == 0:
                self.logger.error(
                    '{} has an invalid semi-major axis'.format(desg))
                continue

            tab = tab[:1]  # in case of multiple orbits

            e = tab['Epoch'][0]
            epoch = {'I': '18', 'J': '19', 'K': '20'}[e[0]] + e[1:3] + '-'
            if e[3].isdigit():
                epoch += '0' + e[3] + '-'
            else:
                epoch += str(ord(e[3]) - 55) + '-'
            if e[4].isdigit():
                epoch += '0' + e[4]
            else:
                epoch += str(ord(e[4]) - 55)

            orbit = Orbit.from_dict({
                'targetname': [desg],
                'a': tab['a'].data * u.au,
                'e': tab['e'].data,
                'i': tab['Incl.'].data * u.deg,
                'Omega': tab['Node'].data * u.deg,
                'w': tab['Peri.'].data * u.deg,
                'M': tab['M'].data * u.deg,
                'epoch': [Time(epoch).jd],
                'timescale': ['UTC'],
                'H': tab['H'].data,
                'G': tab['G'].data
            })

            eph = Ephem.from_oo(orbit, epochs=epochs)
            eph.add_column(
                Column(Time(eph['MJD'], format='mjd').jd * u.day), name='Date')
            self.logger.info(desg)
            n, found, _tab = self.find_by_ephemeris(eph)
            summary.append(_tab)

        if len(summary) == 0:
            return self.observation_summary([])
        else:
            return vstack(summary)

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

# 0         1         2         3         4         5         6         7         8         9         10        11        12        13        14        15
# 0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
# Object   H     G    Epoch    M         Peri.      Node       Incl.        e           n         a                     NObs NOpp   Arc    r.m.s.       Orbit ID
# ZVA1B24 15.0  0.15  K18AT 359.86374  359.12333   78.49928   68.71450  0.9648993  0.00113082  91.2446974                101   1   25 days 0.59         NEOCPNomin
