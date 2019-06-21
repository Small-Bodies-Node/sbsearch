# Licensed with the 3-clause BSD license.  See LICENSE for details.
import sqlite3
from itertools import repeat, groupby
from logging import ERROR

import numpy as np
import scipy.ndimage as nd
from astropy.io import ascii
from astropy.time import Time
from astropy.table import Table, Column, vstack
import astropy.units as u
from sbpy.data import Orbit, Ephem

from . import logging, util, schema, ephem
from .schema import Obj, Eph, Obs, Found
from .util import RADec, FieldOfView, Line, Point
from .db import SBDB
from .config import Config


class SBSearch:
    """Search for small Solar System objects in survey data.

    Parameters
    ----------
    config : sbsearch.config.Config, optional
        Use this configuration set.

    save_log : bool, optional
        Set to ``True`` to write the log to the log file.

    disable_log : bool, optional
        Set to ``True`` to disable normal logging; also sets
        ``save_log=False``.

    test : bool, optional
        Create and connect to a test database (ignores
        `config['database']`).

    **kwargs
        If ``config`` is ``None``, pass these additional keyword
        arguments to ``Config`` initialization.

    """

    VMAX = 25  # default magnitude limit for searches

    def __init__(self, config=None, save_log=False, disable_log=False,
                 test=False, **kwargs):
        self.config = Config(**kwargs) if config is None else config
        self.config.update(**kwargs)

        if disable_log:
            save_log = False
            level = ERROR
        else:
            level = None

        fn = self.config['log'] if save_log else '/dev/null'
        self.logger = logging.setup(filename=fn, level=level)

        if test:
            self.db = SBDB.create_test_db()
        else:
            self.db = SBDB(self.config['database'])

        self.verify_database()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.logger.info('Closing database.')
        self.db.session.commit()
        self.db.session.close()
        self.logger.info(Time.now().iso + 'Z')

    def add_found(self, *args, **kwargs):
        location = kwargs.pop('location', self.config['location'])
        args = args + (location,)
        return self.db.add_found(*args, **kwargs)
    add_found.__doc__ = SBDB.add_found.__doc__
    add_found.__doc__.replace('location : string',
                              'location : string, optional')

    def add_found_by_id(self, *args, **kwargs):
        location = kwargs.pop('location', self.config['location'])
        args = args + (location,)
        return self.db.add_found_by_id(*args, **kwargs)
    add_found_by_id.__doc__ = SBDB.add_found_by_id.__doc__
    add_found_by_id.__doc__.replace('location : string',
                                    'location : string, optional')

    def add_observations(self, observations, update=False):
        n = self.db.add_observations(observations, update=update)

        if n < len(observations):
            loglevel = self.logger.warning
        else:
            loglevel = self.logger.info

        if update:
            loglevel('Added or updated {} of {} observations.'
                     .format(n, len(observations)))
        else:
            loglevel('Added {} of {} observations.'
                     .format(n, len(observations)))

    add_observations.__doc__ = SBDB.add_observations.__doc__

    def clean_ephemeris(self, objects, start=None, stop=None):
        """Remove ephemeris between dates (inclusive).

        Parameters
        ----------
        objects : list
            Object IDs or designations.

        start, stop : float or `~astropy.time.Time`, optional
            Date range (inclusive), or ``None`` for unbounded.

        """

        total = 0

        jd_start, jd_stop = util.epochs_to_jd([start, stop])

        self.logger.debug('Removing ephemeris rows:')
        for objid, desg in self.db.resolve_objects(objects):
            if objid is None:
                self.logger.warning('{} not in object table'.format(desg))
                continue

            n = self.db.clean_ephemeris(objid, jd_start, jd_stop)
            self.logger.debug('* {}: {}'.format(desg, n))
            total += n
        self.logger.info('{} ephemeris rows removed.'.format(total))

    def clean_found(self, objects=None, start=None, stop=None):
        """Remove found entries between dates (inclusive).

        Parameters
        ----------
        objects : list, optional
            Object IDs or designations, or ``None`` to remove all
            found entries.

        start, stop : float or `~astropy.time.Time`, optional
            Date range (inclusive), or ``None`` for unbounded.

        Returns
        -------
        count : int
            Number of rows removed.

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
        return total

    def find_by_ephemeris(self, eph):
        """Find object based on ephemeris.

        Use for objects that are not in the database.

        Parameters
        ----------
        eph : `~sbpy.data.Ephem`
            Ephemeris of object to find, requires RA, Dec, and date
            columns.

        Returns
        -------
        obsids : list
            Observation IDs with this object.

        """

        if len(eph) == 1:
            raise ValueError('Cannot search single-point ephemerides.')

        # make sure ephemeris is in time order
        ordered_eph = Ephem.from_table(eph[np.argsort(eph['date'].value)])
        coords = RADec(ordered_eph['RA'], ordered_eph['Dec'])
        jd = ordered_eph['Date'].value

        # search for each ephemeris segment
        found = []
        for i in range(len(ordered_eph) - 1):
            e = Ephem.from_table(ordered_eph[i:i+2])
            target = str(Line.from_ephem(e))

            for obs in self.db.get_observations_intersecting(
                    target, start=jd[i], stop=jd[i+1]).all():
                # interpolate to observation time
                c = util.spherical_interpolation(
                    coords[i], coords[i+1], jd[i], jd[i+1],
                    (obs.jd_start + obs.jd_stop) / 2)

                # second pass
                if self.db.observation_covers(obs, c):
                    found.append(obs)

        self.logger.info('{} observations found.'.format(len(found)))
        obsids = [f.obsid for f in found]

        return obsids

    def find_object(self, obj, start=None, stop=None, vmax=VMAX,
                    source=None, location=None, save=True, update=False,
                    **kwargs):
        """Find observations covering object.

        Parameters
        ----------
        obj : int or str
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

        location : string, optional
            Observer location.  If ``None``, use config default.

        save : bool, optional
            Save observations to found database.

        update : bool, optional
            If ``True``, update metadata for already found objects.

        **kwargs
            Additional keyword arguments for `~ephem.generate`.

        Returns
        -------
        obsids : list
            Observation IDs with this object.

        foundids : list
            All found IDs.

        newids : list
            Only the newly found IDs.

        """

        objid, desg = self.db.resolve_object(obj)
        if objid is None:
            objid = self.db.add_object(desg)

        obs_range = self.db.get_observation_date_range()
        if start is None:
            start = obs_range[0]
        if stop is None:
            stop = obs_range[-1]
        t_start, t_stop = util.epochs_to_time((start, stop))
        jd_start, jd_stop = util.epochs_to_jd((start, stop))

        location = self.config['location'] if location is None else location

        if source is not None:
            epochs = {
                'start': t_start.iso[:16],
                'stop': t_stop.iso[:16],
                'step': None
            }
            eph = ephem.generate(desg, location, epochs, source=source,
                                 **kwargs)
            obsids = []

            # group ephemerides into segments based on V magnitude, and
            # search the database for each segment
            eph['_v_'] = util.vmag_from_eph(eph)
            groups = groupby(eph, lambda e: e['_v_'] <= vmax)
            for bright_enough, group in groups:
                if bright_enough:
                    group_eph = Ephem.from_table(vstack(list(group)))
                    obsids.extend(self.find_by_ephemeris(group_eph))
        else:
            eph = (self.db.get_ephemeris(objid, jd_start, jd_stop)
                   .filter(Eph.vmag <= vmax)
                   .all())
            obsids = []
            if len(eph) > 0:
                target = str(util.Line.from_eph(eph))
                obs = self.db.get_observations_intersecting(
                    target, start=jd_start, stop=jd_stop).all()
                obsids = [o.obsid for o in obs]

        foundids = []
        newids = []
        if save and len(obsids) > 0:
            foundids, newids = self.db.add_found_by_id(
                objid, obsids, location, update=update,
                cache=kwargs.get('cache', False))

        return obsids, foundids, newids

    def find_objects(self, objects, start=None, stop=None, progress=None,
                     **kwargs):
        """Find observations covering any of these objects.

        Parameters
        ----------
        objects : list or ``None``
            Objects to find, designations and/or object IDs.  ``None``
            to search for all objects.

        start : float or `~astropy.time.Time`, optional
            Search after this epoch, inclusive.

        stop : float or `~astropy.time.Time`, optional
            Search before this epoch, inclusive.

        progress : ProgressWidget, optional
            Report progress to this `~logging` widget.

        **kwargs
            Keyword arguments for `~find_object`.

        Returns
        -------
        obsids : list
            Observation IDs with this object.

        foundids : list
            All found IDs.

        newids : list
            Only the newly found IDs.

        """

        _objects = self.db.resolve_objects(objects)
        jd_start, jd_stop = util.epochs_to_jd((start, stop))

        if start is None and stop is None:
            s = 'in all observations'
        elif start is None:
            t = util.epochs_to_time([stop])[0]
            s = 'in observations ending {} UT'.format(t.iso[:16])
        elif stop is None:
            t = util.epochs_to_time([start])[0]
            s = 'in observations starting {} UT'.format(t.iso[:16])
        elif start == stop:
            t = util.epochs_to_time([start])[0]
            s = 'in observations on {}'.format(t.iso[:10])
        else:
            t = util.epochs_to_time([start, stop])
            s = 'in observations from {} UT to {} UT'.format(
                t[0].iso[:16], t[1].iso[:16])

        s += ', V<={:.1f}'.format(kwargs.get('vmax', self.VMAX))

        self.logger.info('Searching for {} object{} {}.'.format(
            len(_objects), '' if len(_objects) == 1 else 's', s))

        obsids = []
        foundids = []
        newids = []
        progress = logging.ProgressTriangle(1, self.logger, base=2)
        for objid, desg in _objects:
            o, f, n = self.find_object(objid, start=jd_start, stop=jd_stop,
                                       **kwargs)
            obsids.extend(o)
            foundids.extend(f)
            newids.extend(n)

            self.logger.debug('* {} x{}, {} saved'.format(
                desg, len(obsids), len(n)))
            progress.update(len(obsids))

        progress.done()
        self.logger.info(
            'Found in {} observations ({} saved).'.format(
                len(obsids), len(newids)))

        return obsids, foundids, newids

    def list_objects(self):
        """List available objects.

        Returns
        -------
        tab : `~astropy.table.Table`

        """
        objects = [(obj.objid, obj.desg) for obj in
                   self.db.get_objects().all()]
        tab = Table(rows=objects,
                    names=('object ID', 'designation'))
        return tab

    def summarize_found(self, objects=None, start=None, stop=None):
        """Summarize found objects.

        Parameters
        ----------
        objects : list, optional
            Object IDs or desginations to analyze.

        start, stop : string or `~astropy.time.Time`, optional
            Date range to search, inclusive.  ``None`` for unbounded
            limit.

        Returns
        -------
        tab : `~astropy.table.Table` or ``None``

        """

        found = self.db.session.query(Found)
        found = util.filter_by_date_range(found, start, stop, Found.jd)
        if objects is not None:
            objids = [obj[0] for obj in self.db.resolve_objects(objects)]
            found = found.filter(Found.objid.in_(objids))

        rows = []
        cols = ('foundid', 'objid', 'obsid', 'jd', 'ra', 'dec', 'dra',
                'ddec', 'unc_a', 'unc_b', 'unc_theta', 'vmag', 'rh',
                'rdot', 'delta', 'phase', 'selong', 'sangle', 'vangle',
                'trueanomaly', 'tmtp')
        for row in found:
            rows.append([getattr(row, k) for k in cols])

        if len(rows) == 0:
            return None

        tab = Table(rows=rows, names=cols)

        return tab

    def summarize_object_coverage(self, cov_type, objects=None, start=None,
                                  stop=None, length=60, source=None):
        """Ephemeris or found coverage for objects.

        Parameters
        ----------
        cov_type: string
            'eph' for ephemeris coverage, 'found' for found coverage.

        objects : list, optional
            Object IDs or desginations to analyze.

        start, stop: string or astropy.time.Time, optional
            Start and stop epochs, parseable by `~astropy.time.Time`.
            Defaults for ``cov_type='eph'`` are min/max ephemeris
            dates, and for ``'found'`` the min/max observation dates.

        length : int, optional
            Length of the strings to create.

        source : string, optional
            Limit found observations to this obs table source.

        Returns
        -------
        tab : `~astropy.table.Table` or ``None``

        """

        if cov_type not in ['eph', 'found']:
            raise ValueError(
                'coverage type must be eph or found: {}'.format(
                    cov_type))

        if objects is None:
            objects = [(obj.objid, obj.desg)
                       for obj in self.db.get_objects().all()]
        else:
            objects = self.db.resolve_objects(objects)

        jd_start, jd_stop = util.epochs_to_jd((start, stop))

        if None in [jd_start, jd_stop]:
            objids = list([obj[0] for obj in objects])
            if cov_type == 'eph':
                jd_range = self.db.get_ephemeris_date_range(objids=objids)
            elif cov_type == 'found':
                jd_range = self.db.get_observation_date_range(source=source)

            if jd_start is None:
                jd_start = jd_range[0]

            if jd_stop is None:
                jd_stop = jd_range[1]

        bins = np.linspace(jd_start, jd_stop, length + 1)

        coverage = []
        for objid, desg in objects:
            if objid is None:
                coverage.append(('', desg, '-' * length))
                continue

            if cov_type == 'eph':
                query = self.db.get_ephemeris(objid, jd_start, jd_stop)
                jd = np.array([q.jd for q in query])
            elif cov_type == 'found':
                query = self.db.get_found(obj=objid, start=jd_start,
                                          stop=jd_stop)
                jd = np.array([q.jd for q in query])

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

    def summarize_observations(self, obsids):
        """Summarize observations.

        Parameters
        ----------
        obsids : tuple or list of int
            Observation IDs to summarize.

        Returns
        -------
        summary : `~astropy.table.Table` or ``None``
            Summary table.

        """

        obs = self.db.session.query(
            Obs.obsid, Obs.source, Obs.jd_start, Obs.jd_stop,
            Obs.fov.ST_AsGeoJSON(), Obs.seeing, Obs.airmass,
            Obs.maglimit).filter(Obs.obsid.in_(obsids)).all()
        names = ('obsid', 'source', 'start', 'stop',
                 'FOV', 'seeing', 'airmass', 'maglimit')

        if len(obs) == 0:
            return None
        else:
            return Table(rows=obs, names=names)

    def update_ephemeris(self, objects, start, stop, step=None,
                         source='jpl', cache=False):
        """Update object ephemeris table.

        Parameters
        ----------
        objects: list of string or integer
            Designations or object IDs.

        start, stop: string or astropy.time.Time
            Start and stop epochs, parseable by `~astropy.time.Time`.
            Previously defined ephemerides will be removed.

        step: string or astropy.unit.Quantity, optional
            Integer step size in hours or days.  If `None`, then an
            adaptable step size will be used, based on distance to the
            observer.

        source : string
            Source to use: 'mpc' or 'jpl'.

        cache : bool, optional
            Use cached ephemerides; primarily for testing.

        """

        jd_start, jd_stop = util.epochs_to_jd([start, stop])

        cleaned = 0
        added = 0
        for obj in objects:
            objid, desg = self.db.resolve_object(obj)
            self.logger.debug('* {}'.format(desg))
            if objid is None:
                objid = self.db.add_object(desg)

            n = self.db.clean_ephemeris(objid, jd_start, jd_stop)
            cleaned += n

            n = self.db.add_ephemeris(
                objid, self.config['location'], jd_start, jd_stop, step=step,
                source=source, cache=cache)
            added += n

        self.logger.info('{} rows deleted from eph table'.format(cleaned))
        self.logger.info('{} rows added to eph table'.format(added))

    def verify_database(self, names=[]):
        """Verify database."""
        self.db.verify_database(self.logger, names=names)

# 0         1         2         3         4         5         6         7         8         9         10        11        12        13        14        15
# 0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
# Object   H     G    Epoch    M         Peri.      Node       Incl.        e           n         a                     NObs NOpp   Arc    r.m.s.       Orbit ID
# ZVA1B24 15.0  0.15  K18AT 359.86374  359.12333   78.49928   68.71450  0.9648993  0.00113082  91.2446974                101   1   25 days 0.59         NEOCPNomin
