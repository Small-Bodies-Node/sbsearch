# Licensed with the 3-clause BSD license.  See LICENSE for details.
from itertools import groupby
from logging import ERROR, DEBUG

import numpy as np
from astropy.time import Time
from astropy.table import Table
from geoalchemy2.functions import ST_MakeLine

from . import logging, util, ephem
from .schema import Eph, Obs, Found
from .util import RADec, Line
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

    session : sqlalchemy Session, optional
        Ignore ``config['database']`` and use this session instead.

    **kwargs
        If ``config`` is ``None``, pass these additional keyword
        arguments to ``Config`` initialization.

    """

    VMAX = 25  # default magnitude limit for searches

    def __init__(self, config=None, logger_name='SBSearch', save_log=False,
                 disable_log=False, test=False, session=None, debug=False, **kwargs):
        self.config = Config(**kwargs) if config is None else config
        self.config.update(**kwargs)
        self.debug = debug

        if disable_log:
            save_log = False
            level = ERROR
        elif self.debug:
            level = DEBUG
        else:
            level = None

        fn = self.config['log'] if save_log else '/dev/null'
        self.logger = logging.setup(filename=fn, level=level, name=logger_name)

        if test:
            self.db = SBDB.create_test_db()
        elif session is None:
            self.db = SBDB(self.config['database'])
        else:
            self.db = SBDB(session)

        self.verify_database()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.db.session.commit()
        self.db.session.close()
        self.logger.debug('Disconnected from database.')

    def add_found(self, *args, **kwargs):
        return self.db.add_found(*args, **kwargs)
    add_found.__doc__ = SBDB.add_found.__doc__

    def add_found_by_id(self, *args, **kwargs):
        return self.db.add_found_by_id(*args, **kwargs)
    add_found_by_id.__doc__ = SBDB.add_found_by_id.__doc__

    def add_observations(self, observations, update=False):
        """Add observations to database.

        If observations already exist for a given observation ID, the
        old data are updated.


        Parameters
        ----------
        observations : list of Obs
            Observations to insert.

        update : bool, optional
            Update database in case of duplicates

        Returns
        -------
        n : int
            Number of inserted or updated rows.

        """
        n = self.db.add_observations(
            observations, update=update, logger=self.logger)

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

    def find_by_ephemeris(self, eph, source=Obs, vmax=VMAX):
        """Find object based on ephemeris.

        Use for objects that are not in the database.

        Parameters
        ----------
        eph : `~sbpy.data.Ephem`
            Ephemeris of object to find, requires RA, Dec, and Date
            columns.  May be an ``Ephem`` object or a list of
            dictionaries.  Must be in time order.

        source : sqlalchemy mapping, optional
            Source for observations.  Use ``schema.Obs`` to search all
            sources.  Otherwise, pass the specific source object.

        vmax : float, optional
            Require epochs brighter than this limit.

        Returns
        -------
        obsids : list
            Observation IDs with this object.

        """

        if len(eph) == 1:
            raise ValueError('Cannot search single-point ephemerides.')

        coords = RADec(np.squeeze(eph['RA']), np.squeeze(eph['DEC']))
        jd = np.squeeze(np.array(eph['jd']))
        v = util.vmag_from_eph(eph)

        self.logger.debug('Initialized')

        # Search source observations for all ephemeris segments
        n = len(coords)
        found = []
        chunk = 10  # limit the intersection query to this many segments at a time
        for i in range(0, n - 1, chunk):
            j = min(i + chunk, n - 1)
            segments = ST_MakeLine([
                str(Line(coords[k:k+2])) for k in range(i, j)
            ])
            _jd = jd[i:j+1]
            observations = self.db.get_observations_intersecting(
                segments, start=_jd[0], stop=_jd[-1], source=source
            ).all()

            # now we have observations that interset the target's
            # ephemeris; second pass is to check the precise position
            for obs in observations:
                # interpolate to observation time
                jd_mid = (obs.jd_start + obs.jd_stop) / 2
                k = np.searchsorted(jd, jd_mid, side='left')

                if not np.any(v[k-1:k+1] < vmax):
                    continue

                c = util.spherical_interpolation(
                    coords[k-1], coords[k], jd[k-1], jd[k],
                    jd_mid)

                if self.db.observation_covers(obs, c):
                    found.append(obs)

        # eliminate duplicates
        found = list(set(found))

        n = len(found)
        self.logger.info('Found {} observation{} in {}'.format(
            n, '' if n == 1 else 's', source.__data_source_name__))

        obsids = [f.obsid for f in found]
        return obsids

    def find_object(self, obj, start=None, stop=None, source=Obs, vmax=VMAX,
                    eph_source=None, save=True, update=False, **kwargs):
        """Find observations covering object.

        Parameters
        ----------
        obj : int or str
            Object to find.

        start : float or `~astropy.time.Time`, optional
            Search after this epoch, inclusive.

        stop : float or `~astropy.time.Time`, optional
            Search before this epoch, inclusive.

        source : sqlalchemy mapping, optional
            Source for observations.  Use ``schema.Obs`` to search all
            sources.  Otherwise, pass the specific source object.

        vmax : float, optional
            Require epochs brighter than this limit.

        eph_source : string, optional
            Ephemeris source: ``None`` for internal database, 'mpc' or
            'jpl' for online ephemeris generation.

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
        self.logger.info('Searching for {} in {}'.format(
            desg, source.__data_source_name__))

        obs_range = self.db.get_observation_date_range(source=source)
        if start is None:
            start = obs_range[0]
        if stop is None:
            stop = obs_range[-1]
        t_start, t_stop = util.epochs_to_time((start, stop))
        jd_start, jd_stop = util.epochs_to_jd((start, stop))

        self.logger.debug('Set up target and time span')
        if eph_source is not None:
            epochs = {
                'start': t_start.iso[:16],
                'stop': t_stop.iso[:16],
                'step': None
            }
            eph = ephem.generate(desg, source.__obscode__, epochs, source=eph_source,
                                 **kwargs)
            self.logger.debug(
                'Obtained ephemeris from {}'.format(eph_source))

            obsids = self.find_by_ephemeris(eph, source=source, vmax=vmax)
        else:
            eph = (self.db.get_ephemeris(objid, jd_start, jd_stop)
                   .filter(Eph.vmag <= vmax)
                   .all())
            self.logger.debug(
                'Obtained ephemeris from internal database')
            obsids = []
            n = len(eph)
            if n > 0:
                self.logger.warning('Untested search by ephemeris')
                chunk = 10  # limit the intersection query to this many segments at a time
                for i in range(0, n - 1, chunk):
                    j = min(i + chunk, n - 1)
                    segments = ST_MakeLine([
                        str(util.Line.from_eph(eph[k:k+2])) for k in range(i, j)
                    ])
                    obs = self.db.get_observations_intersecting(
                        segments, start=eph[0].jd, stop=eph[-1].jd, source=source
                    ).all()

                    obsids.extend([o.obsid for o in obs])

        self.logger.debug('Completed search')

        foundids = []
        newids = []
        n = len(obsids)
        if save and n > 0:
            foundids, newids = self.db.add_found_by_id(
                objid, obsids, update=update, cache=kwargs.get('cache', False))

        self.logger.debug('Added found observations to database.')
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

        source : sqlalchemy mapping, optional
            Source for observations.  Use ``schema.Obs`` to search all
            sources.  Otherwise, pass the specific source object.

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
                objid, Obs.__obscode__, jd_start, jd_stop, step=step,
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
