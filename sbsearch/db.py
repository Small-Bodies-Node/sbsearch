# Licensed with the 3-clause BSD license.  See LICENSE for details.
import re
import sqlite3
import itertools
from logging import Logger
import struct

import numpy as np
from numpy import pi
from astropy.coordinates import Angle
import astropy.units as u
from astropy.time import Time
from astropy.table import vstack
from sbpy.data import Ephem, Names, Orbit

from . import util, schema
from .util import RADec
from .logging import ProgressTriangle
from .exceptions import BadObjectID, NoEphemerisError, SourceNotFoundError


class SBDB(sqlite3.Connection):
    """Database object for SBSearch."""

    DB_NAMES = ['obj', 'eph', 'eph_tree', 'delete_eph_from_tree',
                'delete_object_from_eph', 'obs', 'obs_sources',
                'obs_tree', 'delete_obs_from_obs_tree', 'found',
                'delete_obs_from_found', 'delete_object_from_found']

    def __init__(self, *args):
        super().__init__(*args)
        sqlite3.register_adapter(np.int64, int)
        sqlite3.register_adapter(np.float64, float)
        self.row_factory = sqlite3.Row

        self.execute('PRAGMA foreign_keys = 1')
        self.execute('PRAGMA recursive_triggers = 1')

    @classmethod
    def make_test_db(cls, target=':memory:'):
        from itertools import repeat
        from .test.skytiles import sky_tiles, N_tiles

        db = sqlite3.connect(target, 5, 0, "DEFERRED", True, SBDB)
        db.verify_database(Logger('test'))
        db.add_object('C/1995 O1')
        db.add_object('2P')

        obsids = range(N_tiles**2)
        start = 2458119.5 + np.arange(N_tiles**2) * 30 / 86400
        stop = start + 30 / 86400
        columns = [obsids, repeat('test'), start, stop, sky_tiles]
        db.add_observations(zip(*columns))
        db.add_ephemeris(2, '500', 2458119.5, 2458121.5, step='1d',
                         source='jpl', cache=True)
        db.add_found_by_id(2, [1, 2, 3], '500', cache=True)
        return db

    def verify_database(self, logger, names=[], script=''):
        """Verify SBSearch tables, triggers, etc.

        Parameters
        ----------
        logger : `~logging.Logger`
            Log messages to this logger.

        names : list, optional
            Additional database names to consider.

        script : string, optional
            Additional SQL commands to execute if any names are missing.

        """

        c = self.execute("SELECT name FROM sqlite_master")
        existing_names = list([row[0] for row in c])
        missing = False
        for name in self.DB_NAMES + names:
            if name not in existing_names:
                missing = True
                logger.error('{} is missing from database'.format(name))

        if missing:
            schema.create(self)
            if len(script) > 0:
                self.executescript(script)
            logger.info('Created database tables and triggers.')

    def add_ephemeris(self, objid, location, jd_start, jd_stop, step=None,
                      source='jpl', cache=False):
        """Add ephemeris data to databse.

        Parameters
        ----------
        objid : int
            Database object ID.

        location : string
            Observer location.

        jd_start, jd_stop : float
            Julian date range (inclusive).

        step : string or `~astropy.units.Quantity`, optional
            Ephemeris step size or ``None`` to use an adaptable size
            based on heliocentric distance.

        source : string
            Source to use: 'mpc' or 'jpl'.

        cache : bool, optional
            Use cached ephemerides; primarily for testing.

        Returns
        -------
        count : int
            Number of inserted rows.

        """

        if step is None:
            # Adaptable time step
            # daily ephemeris for delta > 1
            count = self.add_ephemeris(objid, location, jd_start, jd_stop,
                                       step='1d', source=source,
                                       cache=cache)
            # ZChecker analysis, Oct 2018: error ~ 2" / delta for 6h time step
            for limit, substep in ((1, '4h'), (0.25, '1h')):
                eph = self.get_ephemeris(
                    objid, jd_start, jd_stop, columns='jd,delta')
                groups = itertools.groupby(eph, lambda e: e[1] < limit)
                for inside, epochs in groups:
                    if not inside:
                        continue

                    jd = list([e[0] for e in epochs])
                    if len(jd) > 1:
                        count -= self.clean_ephemeris(objid, jd[0], jd[-1])
                        count += self.add_ephemeris(
                            objid, location, jd[0], jd[-1], step=substep,
                            source=source, cache=cache)
        else:
            desg = self.resolve_object(objid)[1]
            step = u.Quantity(step)
            count = 0
            total = int(round((jd_stop - jd_start) / step.to('day').value) + 1)
            next_step = jd_start
            today = Time.now().iso[:10]
            while count < total:
                n = total - count

                if source == 'mpc':
                    epochs = {'start': next_step, 'step': step, 'number': n}
                elif source == 'jpl':
                    stop = next_step + step.to('day').value * (n - 1)
                    epochs = {
                        'start': str(Time(next_step, format='jd').iso),
                        'stop': str(Time(stop, format='jd').iso),
                        'step': str(n - 1)
                    }

                eph = self.get_ephemeris_exact(
                    desg, location, epochs, source=source, cache=cache)

                jd = util.epochs_to_jd(eph['Date'].value)
                coords = RADec(eph['RA'], eph['Dec'])
                half_step = step.to('d').value / 2

                vmag = util.vmag_from_eph(eph)

                for i in range(len(eph)):
                    # save to ephemeris table
                    cursor = self.execute('''
                    INSERT OR REPLACE INTO eph
                    VALUES (null,?,?,?,?,?,?,?,?,?,?)
                    ''', (objid,
                          jd[i],
                          eph['r'][i].value,
                          eph['Delta'][i].value,
                          eph['RA'][i].to('rad').value,
                          eph['Dec'][i].to('rad').value,
                          eph['dRA cos(Dec)'][i].to('arcsec/hr').value,
                          eph['ddec'][i].to('arcsec/hr').value,
                          vmag[i],
                          today))

                    # save to ephemeris tree
                    if i == 0:
                        indices = (1, 0, 1)
                    elif i == len(eph) - 1:
                        indices = (-2, -1, -2)
                    else:
                        indices = (i - 1, i, i + 1)

                    c = tuple((coords[j] for j in indices))
                    _jd = tuple((jd[j] for j in indices))
                    # jd to mjd conversion in eph_to_limits
                    limits = util.eph_to_limits(c, _jd, half_step)
                    self.execute('''
                    INSERT OR REPLACE INTO eph_tree VALUES (
                      last_insert_rowid(),?,?,?,?,?,?,?,?
                    )''', limits)

                count += len(eph)
                next_step = jd_start + (step * (count + 1)).to(u.day).value

        return count

    def add_found(self, objid, observations, location, update=False,
                  cache=False):
        """Add found objects to found database.

        Parameters
        ----------
        objid : int
            Found object ID.

        observations : tuple or list
            Observations with found objects.

        location : string
            Observer location.

        update : bool, optional
            ``True`` to overwrite existing entries.

        cache : bool, optional
            Use cached ephemerides; primarily for testing.

        Returns
        -------
        foundids : list
            New found IDs.  If ``update`` is ``False``, found IDs that
            already exist will not be returned.

        """

        observations = list(observations)
        foundids = np.zeros(len(observations))

        # already in found database?
        if not update:
            missing = []
            obs_missing = []
            for i, obs in enumerate(observations):
                f = self.execute('''
                SELECT foundid FROM found
                WHERE objid=? AND obsid=?
                ''', [objid, obs['obsid']]).fetchone()
                if f:
                    foundids[i] = f[0]
                else:
                    missing.append(i)
                    obs_missing.append(obs)

            if len(missing) > 0:
                foundids[missing] = self.add_found(
                    objid, obs_missing, location, update=True,
                    cache=cache)

            return foundids[missing]

        jd = np.array([(obs['jd_start'] + obs['jd_stop']) / 2
                       for obs in observations])

        jd_sorted, unsort_jd = np.unique(jd, return_inverse=True)
        eph = self.get_ephemeris_exact(objid, location, jd_sorted,
                                       source='jpl', cache=cache)
        orb = self.get_orbit_exact(objid, jd_sorted, cache=cache)

        # restore requested order
        eph = Ephem.from_table(eph[unsort_jd])
        orb = Ephem.from_table(orb[unsort_jd])

        vmag = util.vmag_from_eph(eph)
        sangle = Angle(eph['sunTargetPA'] - 180 * u.deg)
        sangle = sangle.wrap_at(360 * u.deg).deg
        vangle = Angle(eph['velocityPA'] - 180 * u.deg)
        vangle = vangle.wrap_at(360 * u.deg).deg
        Tp = Time(orb['Tp_jd'], format='jd', scale='tt').utc.jd
        tmtp = Tp - jd

        for i, obs in enumerate(observations):
            data = [objid, obs['obsid'], jd[i]]
            for col, unit in (('ra', 'deg'), ('dec', 'deg'),
                              ('dRA cos(Dec)', 'arcsec/hr'),
                              ('ddec', 'arcsec/hr'),
                              ('RA_3sigma', 'arcsec'),
                              ('DEC_3sigma', 'arcsec')):
                data.append(eph[col][i].to(unit).value)
            data.append(vmag[i])
            for col, unit in (('r', 'au'), ('r_rate', 'km/s'),
                              ('delta', 'au'), ('alpha', 'deg'),
                              ('elong', 'deg')):
                data.append(eph[col][i].value)
            data.extend((sangle[i], vangle[i], orb['nu'][i].value,
                         tmtp[i]))

            c = self.execute('''
            INSERT OR REPLACE INTO found
            VALUES (null,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
            ''', data)
            foundids[i] = c.lastrowid

        self.commit()

        return foundids

    def add_found_by_id(self, objid, obsids, location, **kwargs):
        """Add found objects to found database using observation ID.

        Parameters
        ----------
        objid : int
            Found object ID.

        obsids : tuple or list of int
            Observation IDs with found objects.

        location : string
            Observer location.

        **kwargs
            Any ``add_found`` keyword arguments.

        """

        if isinstance(obsids[0], sqlite3.Row):
            _obsids = [row[0] for row in obsids]
        else:
            _obsids = obsids
        observations = self.get_observations_by_id(_obsids)
        return self.add_found(objid, observations, location, **kwargs)

    def add_object(self, desg):
        """Add new object to object database.

        Parameters
        ----------
        desg : string
            Object designation.

        Returns
        -------
        objid : int
            The new object ID.

        """

        if not isinstance(desg, str):
            raise ValueError('desg must be a string')

        c = self.execute('''INSERT INTO obj VALUES (null,?)''', [desg])
        return c.lastrowid

    def add_observations(self, rows, other_cmd=None, other_rows=None,
                         logger=None):
        """Add observations to observation table.

        RA, Dec must be in radians.  If observations already exist for
        a given observation ID, the new data are ignored.

        Parameters
        ----------
        rows : iterable
            Rows of data to insert: obsid (or None), source table,
            start JD, stop JD, points, where points is an array of RA
            and Dec points in radians describing the observation field
            of view.  See notes for a detailed description.

        other_cmd : string, optional
            Additional insert command to run after each row, e.g., to
            insert rows into another observation table.  Use
            last_row_id() to keep obsids synced.

        other_rows : iterable
            Parameters for ``other_cmd``.

        logger : `~logging.Logger`, optional
            Report progress to this logger.

        Notes
        -----
        The points array is N * 10 elements long, where N is the
        number of rectangles to describe.  Each rectangle is described
        by 5 points, the center RA and Dec, followed by 4 corners:
        ``[ra_center, dec_center, ra1, dec1, ra2, dec2, ra3, dec3,
        ra4, dec4]``.  The order of the corners is not important.

        """
        def row_iterator(rows, other_rows):
            if other_rows:
                for row, other in zip(rows, other_rows):
                    yield row, other
            else:
                for row in rows:
                    yield row, None

        obs_cmd = 'INSERT OR IGNORE INTO obs VALUES (?,?,?,?,?)'

        tree_cmd = '''
        INSERT OR IGNORE INTO obs_tree VALUES (?,?,?,?,?,?,?,?,?)
        '''

        if logger is not None:
            tri = ProgressTriangle(1, logger, base=10)

        for row, other_row in row_iterator(rows, other_rows):
            fov = struct.pack('10d', *list(row[4]))
            c = self.execute(obs_cmd, (row[0], row[1], row[2], row[3], fov))
            obsid = c.lastrowid

            if other_cmd:
                self.execute(other_cmd, other_row)

            mjd = (np.array((row[2], row[3])) - 2400000.5)
            ra = row[4][::2]
            dec = row[4][1::2]
            xyz = util.rd2xyz(ra, dec)
            self.execute(
                tree_cmd, [obsid, min(mjd), max(mjd),
                           min(xyz[0]), max(xyz[0]),
                           min(xyz[1]), max(xyz[1]),
                           min(xyz[2]), max(xyz[2])])

            if logger is not None:
                tri.update()

    def clean_ephemeris(self, objid, jd_start, jd_stop):
        """Remove ephemeris between dates (inclusive).

        Parameters
        ----------
        objid: int
            Object ID.

        jd_start, jd_stop: float or ``None``
            Julian date range(inclusive), or ``None`` for unbounded.

        Returns
        -------
        n: int
            Number of removed rows.

        """

        constraints = [('objid=?', objid)]
        constraints.extend(util.date_constraints(jd_start, jd_stop))
        cmd, parameters = util.assemble_sql(
            'DELETE FROM eph', [], constraints)

        c = self.execute(cmd, parameters)
        return c.rowcount

    def clean_found(self, objid, jd_start, jd_stop):
        """Remove found objects.

        Parameters
        ----------
        objid: int
            Object ID or ``None`` for all objects.

        jd_start, jd_stop: float or ``None``
            Julian date range(inclusive), or ``None`` for unbounded.

        Returns
        -------
        count: int
            Number of rows removed.

        """

        constraints = []
        if objid:
            constraints.append(('objid=?', objid))

        constraints.extend(
            util.date_constraints(jd_start, jd_stop, 'obsjd'))

        cmd, parameters = util.assemble_sql(
            'DELETE FROM found', [], constraints)

        c = self.execute(cmd, parameters)
        return c.rowcount

    def get_ephemeris(self, objid, jd_start, jd_stop, columns='*',
                      generator=False, order=True):
        """Get ephemeris data from database.

        Parameters
        ----------
        objid : int or ``None``
            Database object ID or ``None`` for any object.

        jd_start, jd_stop : float
            Julian date range (inclusive).

        columns : string, optional
            Columns to return.  Default all.

        generator : bool, optional
            ``True`` to return a generator.

        order : bool, optional
            ``True`` to sort by Julian date.

        Returns
        -------
        eph : generator or list of rows

        """

        cmd = 'SELECT {} FROM eph'.format(columns)

        constraints = []
        parameters = []
        if objid is None:
            constraints.append(('objid NOTNULL', None))
        else:
            constraints.append(('objid=?', objid))

        constraints.extend(util.date_constraints(jd_start, jd_stop))
        cmd, parameters = util.assemble_sql(cmd, parameters, constraints)

        if order:
            cmd += ' ORDER BY jd'

        c = self.execute(cmd, parameters)
        if generator:
            return util.iterate_over(c)
        else:
            return c.fetchall()

    def get_ephemeris_exact(self, obj, location, epochs, source='jpl',
                            orbit=None, cache=False):
        """Generate ephemeris at specific epochs from external source.

        Parameters
        ----------
        obj : int
            Object designation or object ID; not required to be in
            database.

        location : string
            Observer location.

        epochs : array-like or dict
            Compute ephemeris at these epochs.  For arrays, must be
            floats (for Julian date) or else parsable by
            `~astropy.time.Time`.  Dictionaries are passed to the
            ephemeris source as is.

        source : string, optional
            Source to use: 'mpc', 'jpl', or 'oorb'.  'oorb' requires
            ``orbit`` parameter.

        orbit : `~sbpy.data.Orbit`, optional
            Orbital elements for ``source=oorb``.

        cache : bool, optional
            Use cached ephemerides; primarily for testing.

        Returns
        -------
        eph : `~sbpy.data.ephem.Ephem`
            Ephemeris.

        """

        if source not in ['mpc', 'jpl', 'oorb']:
            raise ValueError('Source must be "mpc" or "jpl".')

        if isinstance(obj, str):
            desg = obj
        else:
            desg = self.resolve_object(obj)[1]

        if isinstance(epochs, dict):
            _epochs = epochs
        else:
            _epochs = util.epochs_to_jd(epochs)

            d = np.diff(_epochs)
            if any(d <= 0):
                raise ValueError(
                    'Epoch dates must be increasing and unique: {}'.format(
                        _epochs))

            if len(_epochs) > 300:
                eph = None
                N = np.ceil(len(_epochs) / 200)
                for e in np.array_split(_epochs, N):
                    _eph = self.get_ephemeris_exact(
                        obj, location, e, source=source, cache=cache)
                    if eph:
                        eph.add_rows(_eph)
                    else:
                        eph = _eph
                return eph

        if source == 'mpc':
            eph = Ephem.from_mpc(desg, epochs=_epochs,
                                 location=location,
                                 proper_motion='sky',
                                 proper_motion_unit='rad/s',
                                 cache=cache)
        elif source == 'jpl':
            kwargs = dict(epochs=_epochs,
                          location=location,
                          quantities='1,3,8,9,19,20,23,24,27,36',
                          cache=cache)
            if Names.asteroid_or_comet(desg) == 'comet':
                kwargs.update(id_type='designation',
                              closest_apparition=True,
                              no_fragments=True)

            eph = Ephem.from_horizons(desg, **kwargs)
        elif source == 'oorb':
            eph = Ephem.from_oo(orbit, epochs=_epochs, location=location)

        return eph

    def get_ephemeris_date_range(self, objids=None):
        """Ephemeris date limits.


        Parameters
        ----------
        objids : list, optional
            Limit query to these object IDs.


        Returns
        -------
        jd_min, jd_max : float

        """

        if objids:
            q = ','.join('?' * len(objids))
            cmd = '''
            SELECT MIN(jd) FROM eph INNER JOIN (
              SELECT ephid FROM eph_tree WHERE mjd0 <= (
                SELECT MIN(mjd1) FROM eph_tree
                INNER JOIN eph USING (ephid)
                WHERE objid IN ({})
              )
            ) USING (ephid)
            '''.format(q)
            jd_min = self.execute(cmd, objids).fetchone()[0]

            cmd = '''
            SELECT MAX(jd) FROM eph INNER JOIN (
              SELECT ephid FROM eph_tree WHERE mjd1 >= (
                SELECT MAX(mjd0) FROM eph_tree
                INNER JOIN eph USING (ephid)
                WHERE objid IN ({})
              )
            ) USING (ephid)
            '''.format(q)
            jd_max = self.execute(cmd, objids).fetchone()[0]
        else:
            cmd = '''
            SELECT MIN(jd) FROM eph INNER JOIN (
              SELECT ephid FROM eph_tree WHERE mjd0 <= (
                SELECT MIN(mjd1) FROM eph_tree
              )
            ) USING (ephid)
            '''
            jd_min = self.execute(cmd).fetchone()[0]

            cmd = '''
            SELECT MAX(jd) FROM eph INNER JOIN (
              SELECT ephid FROM eph_tree WHERE mjd1 >= (
                SELECT MAX(mjd0) FROM eph_tree
              )
            ) USING (ephid)
            '''
            jd_max = self.execute(cmd).fetchone()[0]

        if jd_min is None or jd_max is None:
            raise NoEphemerisError('No ephemerides for {}'.format(objids))

        return jd_min, jd_max

    def get_ephemeris_interp(self, objid, epochs):
        """Get ephemeris at specific epochs by interpolation.

        Parameters
        ----------
        objid : int
            Object ID.

        epochs : array-like
            Compute ephemeris at these epochs.  Must be floats (for
            Julian date) or parsable by `~astropy.time.Time`.

        Returns
        -------
        eph : RADec
            Ephemeris.
        vmag : ndarray
            Apparent magnitude estimate.

        """

        ra, dec, vmag = [], [], []
        for jd0 in util.epochs_to_jd(epochs):
            # get two nearest points to epoch
            rows = self.execute('''
            SELECT jd,ra,dec,vmag,ABS(jd - ?) AS dt FROM eph
            WHERE objid=?
              AND dt < 5
            ORDER BY dt LIMIT 2
            ''', [jd0, objid])

            jd, r, d, v, dt = list(zip(*rows))

            i = np.argmin(jd)
            j = np.argmax(jd)

            a = RADec(r[i], d[i], unit='rad')
            b = RADec(r[j], d[j], unit='rad')

            c = util.spherical_interpolation(a, b, jd[i], jd[j], jd0)
            ra.append(c.ra)
            dec.append(c.dec)

            vmag.append(np.interp(jd0, jd, v))

        return RADec(ra, dec, unit='rad'), np.array(vmag)

    def get_ephemeris_segments(self, objid=None, start=None, stop=None,
                               vmax=None):
        """Get ephemeris segments.

        Parameters
        ----------
        objid : int, optional
            Find this object.

        start : float or `~astropy.time.Time`, optional
            Search starting at this epoch, inclusive.

        stop : float or `~astropy.time.Time`, optional
            Search upto and including at this epoch.

        vmax : float, optional
            Require epochs brighter than this limit.

        Returns
        -------
        ephids : ndarray
            Ephemeris IDs for each segment.

        segments : dict of ndarrays
            The segments' x0, x1, y0, etc., suitable for passing to
            ``get_observations_near_box``.

        """

        cmd = '''
        SELECT ephid,mjd0,mjd1,x0,x1,y0,y1,z0,z1 FROM eph_tree
        '''
        constraints = []
        parameters = []

        if objid is not None or vmax is not None:
            cmd += ' INNER JOIN eph USING (ephid)'
            if objid is not None:
                constraints.append(('objid=?', objid))
            if vmax is not None:
                constraints.append(('vmag<=?', vmax))

        if start is not None:
            mjd = util.epochs_to_jd([start])[0] - 2400000.5
            constraints.append(('mjd1 >= ?', mjd))

        if stop is not None:
            mjd = util.epochs_to_jd([stop])[0] - 2400000.5
            constraints.append(('mjd0 <= ?', mjd))

        cmd, parameters = util.assemble_sql(cmd, parameters, constraints)

        keys = ('ephid', 'mjd0', 'mjd1', 'x0', 'x1', 'y0', 'y1', 'z0', 'z1')
        c = self.execute(cmd, parameters)
        values = [np.array(v) for v in zip(*list(c))]

        if len(values) == 0:
            return np.array([]), {}

        segments = dict(zip(keys, values))
        ephids = segments.pop('ephid')
        return ephids, segments

    def get_found(self, obj=None, start=None, stop=None,
                  columns='*', inner_join=None, generator=False):
        """Get found objects by object and/or date range.

        Parameters
        ----------
        obj : int or string, optional
            Find detections of this object.

        start, stop : string or `~astropy.time.Time`, optional
            Date range to search, inclusive.  ``None`` for unbounded
            limit.

        columns : string, optional
            Columns to retrieve, default all.

        inner_join : list, optional
            List of inner join constraints, e.g., `['obs USING
            (obsid)']`.

        generator : bool, optional
            Return a generator rather a list.

        Returns
        -------
        rows : list or generator
            Found object table rows.

        """

        cmd = 'SELECT {} FROM found'.format(columns)
        constraints = util.date_constraints(start, stop, column='obsjd')
        if obj is not None:
            objid = self.resolve_object(obj)[0]
            constraints.append(('objid=?', objid))

        cmd, parameters = util.assemble_sql(cmd, [], constraints,
                                            inner_join=inner_join)
        rows = self.execute(cmd, parameters)

        if generator:
            return util.iterate_over(rows)
        else:
            return list(rows.fetchall())

    def get_found_by_id(self, foundids, columns='*', generator=False):
        """Get found objects by found ID.

        Parameters
        ----------
        foundids : array-like, optional
            Found IDs to retrieve.

        columns : string, optional
            Columns to retrieve, default all.

        generator : bool, optional
            Return a generator rather a list.

        Returns
        -------
        rows : list or generator
            Found object table rows.

        """

        cmd = 'SELECT {} FROM found WHERE foundid=?'.format(columns)

        def g(foundids):
            for foundid in foundids:
                r = self.execute(cmd, [foundid]).fetchone()
                if r is None:
                    return
                else:
                    yield r

        if generator:
            return g(foundids)
        else:
            return list(g(foundids))

    def get_found_by_obsid(self, obsids, columns='*', generator=False):
        """Get found objects by observation ID.

        Parameters
        ----------
        obsids : array-like
            Observation IDs to search.

        columns : string, optional
            Columns to retrieve, default all.

        generator : bool, optional
            Return a generator rather a list.

        Returns
        -------
        rows : list or generator
            Found object table rows.

        """

        cmd = 'SELECT {} FROM found WHERE obsid=?'.format(columns)

        def g(obsids):
            for obsid in obsids:
                r = self.execute(cmd, [obsid]).fetchone()
                if r is None:
                    return
                else:
                    yield r

        if generator:
            return g(obsids)
        else:
            return list(g(obsids))

    def get_objects(self):
        """Return list of all objects.

        Returns
        -------
        objid : ndarray of int
            Object IDs.

        desg : ndarray of string
            Designations.

        """
        objid = []
        desg = []
        for row in util.iterate_over(
                self.execute('SELECT * FROM obj ORDER BY desg+0,desg')):
            objid.append(row[0])
            desg.append(row[1])
        return np.array(objid), np.array(desg)

    def get_observation_date_range(self, source=None):
        """Observation date limits.

        Parameters
        ----------
        survey : string, optional
            Limit query to this observation source.

        Returns
        -------
        jd_min, jd_max : float

        Raises
        ------
        SourceNotFoundError

        """

        if source:
            # minimum
            cmd = '''
            WITH temp AS (
              SELECT obsid,jd_start,mjd0,mjd1 FROM obs_tree
              INNER JOIN obs USING (obsid) WHERE source=?
            )
            SELECT MIN(jd_start) FROM temp WHERE obsid IN (
              SELECT obsid FROM temp WHERE mjd0 <= (
                SELECT MIN(mjd1) FROM temp
              )
            )
            '''
            jd_min = self.execute(cmd, [source]).fetchone()[0]

            # maximum
            cmd = '''
            WITH temp AS (
              SELECT obsid,jd_stop,mjd0,mjd1 FROM obs_tree
              INNER JOIN obs USING (obsid) WHERE source=?
            )
            SELECT MAX(jd_stop) FROM temp WHERE obsid IN (
              SELECT obsid FROM temp WHERE mjd1 >= (
                SELECT MAX(mjd0) FROM temp
              )
            )
            '''
            jd_max = self.execute(cmd, [source]).fetchone()[0]
        else:
            # minimum
            cmd = '''
            SELECT MIN(jd_start) FROM obs INNER JOIN (
              SELECT obsid FROM obs_tree WHERE mjd0 <= (
                SELECT MIN(mjd1) FROM obs_tree
              )
            ) USING (obsid)
            '''
            jd_min = self.execute(cmd).fetchone()[0]

            # maximum
            cmd = '''
            SELECT MAX(jd_stop) FROM obs INNER JOIN (
              SELECT obsid FROM obs_tree WHERE mjd1 >= (
                SELECT MAX(mjd0) FROM obs_tree
              )
            ) USING (obsid)
            '''
            jd_max = self.execute(cmd).fetchone()[0]

        if jd_min is None or jd_max is None:
            raise SourceNotFoundError(
                'No observations for source: {}.'.format(source))

        return jd_min, jd_max

    def get_observations_by_id(self, obsids, columns='*', inner_join=None,
                               generator=False):
        """Get observations by observation ID.

        Parameters
        ----------
        obsids: array-like
            Observation IDs to retrieve.

        columns: string, optional
            Columns to retrieve, default all.

        inner_join : list or tuple of strings, optional
            List of tables and constraints for inner_join, e.g.,
            ``['found USING obsid']``.

        generator: bool, optional
            Return a generator rather than a list.

        Returns
        -------
        rows: list or generator
            Observation table rows.

        """

        _obsids = list(obsids)
        if len(_obsids) == 0:
            return []

        cmd = 'SELECT {} FROM obs'.format(columns)
        q = ','.join(itertools.repeat('?', len(_obsids)))
        constraints = [('obsid IN ({})'.format(q), _obsids)]
        cmd, parameters = util.assemble_sql(cmd, [], constraints,
                                            inner_join=inner_join)
        rows = self.execute(cmd, parameters)

        if generator:
            return util.iterate_over(rows)
        else:
            return list(rows)

    def get_observations_by_date(self, start, stop, columns='*',
                                 inner_join=None, generator=False):
        """Get observations by observation date.

        Parameters
        ----------
        start, stop: string or `~astropy.time.Time`, optional
            Date range to search, inclusive.

        columns: string, optional
            Columns to retrieve, default all.

        inner_join : list or tuple of strings, optional
            List of tables and constraints for inner_join, e.g.,
            ``['found USING obsid']``.

        generator: bool, optional
            Return a generator rather a list.

        Returns
        -------
        rows: list or generator
            Observation table rows.

        """

        jd0, jd1 = util.epochs_to_jd((start, stop))
        constraints = []
        if jd0:
            constraints.append(('mjd1 >= ?', jd0 - 2400000.5))
        if jd1:
            constraints.append(('mjd0 <= ?', jd1 - 2400000.5))

        cmd = '''
        SELECT obsid FROM obs INNER JOIN obs_tree USING (obsid)
        '''
        c = self.execute(*util.assemble_sql(cmd, [], constraints))
        obsids = [row[0] for row in c.fetchall()]

        observations = self.get_observations_by_id(
            obsids, columns=columns, inner_join=inner_join,
            generator=generator)

        return observations

    def get_observations_near(self, ra=None, dec=None, start=None,
                              stop=None, columns='*', inner_join=None,
                              generator=False):
        """Find observations near the given coordinates and/or time.

        Parameters
        ----------
        ra, dec: array-like, float or `~astropy.coordinates.Angle`, optional
            Search this area.  Floats are radians.  Must be at least
            three points.

        start, stop: float or `~astropy.time.Time`, optional
            Search this time range.  Floats are Julian dates.

        columns : string, optional
            Columns to return.

        inner_join : list or tuple of strings, optional
            List of tables and constraints for inner_join, e.g.,
            ``['found USING obsid']``.  The obs table is always
            joined.

        generator: bool, optional
            Return a generator instead of a list.

        Returns
        -------
        obs: list or generator

        """

        if ra is not None:
            if len(ra) < 3:
                raise ValueError('RA requires at least 3 points')

        if dec is not None:
            if len(dec) < 3:
                raise ValueError('Dec requires at least 3 points')

        query = {}
        if ra is not None and dec is None:
            cra = np.cos(ra)
            sra = np.sin(ra)
            query['x0'] = min(cra)
            query['x1'] = max(cra)
            query['y0'] = min(sra)
            query['y1'] = max(sra)
            query['z0'] = -1
            query['z1'] = 1
        elif ra is None and dec is not None:
            sdec = np.sin(dec)
            query['x0'] = -1
            query['x1'] = 1
            query['y0'] = -1
            query['y1'] = 1
            query['z0'] = min(sdec)
            query['z1'] = max(sdec)
        elif ra is not None and dec is not None:
            if len(ra) != len(dec):
                raise ValueError('RA and Dec must have the same length.')
            x, y, z = util.rd2xyz(ra, dec)
            query['x0'] = min(x)
            query['x1'] = max(x)
            query['y0'] = min(y)
            query['y1'] = max(y)
            query['z0'] = min(z)
            query['z1'] = max(z)
        if start is not None:
            mjd = util.epochs_to_time([start]).mjd
            query['mjd0'] = mjd
        if stop is not None:
            mjd = util.epochs_to_time([stop]).mjd
            query['mjd1'] = mjd

        return self.get_observations_near_box(
            inner_join=inner_join, generator=generator, **query)

    def get_observations_near_box(self, columns='*', inner_join=None,
                                  generator=False, **query):
        """Find observations near the given search volume.

        Parameters
        ----------
        mjd0, mjd1, x0, x1, y0, y1, z0, z1 : float or array
            Box(es) to search.

        columns : string, optional
            Columns to return.

        inner_join : list or tuple of strings, optional
            List of tables and constraints for inner_join, e.g.,
            ``['found USING obsid']``.  The obs table is always
            joined.

        generator: bool, optional
            Return a generator instead of a list.

        Returns
        -------
        obs: list or generator

        """

        if len(query) == 0:
            raise ValueError('Nothing to search for.')

        constraints = []
        key2constraint = {
            'mjd0': 'mjd1 >= ?',
            'mjd1': 'mjd0 <= ?',
            'x0': 'x1 >= ?',
            'x1': 'x0 <= ?',
            'y0': 'y1 >= ?',
            'y1': 'y0 <= ?',
            'z0': 'z1 >= ?',
            'z1': 'z0 <= ?'
        }

        # number of boxes to search
        n = max(tuple((np.size(v) for v in query.values())))

        constraints = []
        for k in query.keys():
            constraints.append(key2constraint[k])
        expr = '({})'.format(' AND '.join(constraints))

        parameters = []
        if n == 1:
            for k in query.keys():
                # works for 0- and 1-d arrays, numbers
                parameters.append(float(query[k]))
        else:
            for i in range(n):
                for k in query.keys():
                    parameters.append(query[k][i])

        cmd = 'SELECT {} FROM obs_tree'.format(columns)
        inner_join = [] if inner_join is None else inner_join
        inner_join.append('obs USING (obsid)')
        constraints = [(' OR '.join([expr] * n), parameters)]
        cmd, parameters = util.assemble_sql(cmd, [], constraints,
                                            inner_join=inner_join)
        c = self.execute(cmd, parameters)

        if generator:
            return util.iterate_over(c)
        else:
            return list(c)

    def get_orbit_exact(self, obj, epochs, cache=False):
        """Generate orbital parameters at specific epochs.

        Parameters
        ----------
        obj: string or int
            Object designation or object ID.

        epochs: array-like or dict
            Compute orbital elements at these epochs.  For arrays,
            must be floats(for Julian date) or else parsable by
            `~astropy.time.Time`.

        cache: bool, optional
            Use cached ephemerides; primarily for testing.

        Returns
        -------
        orb: `~sbpy.data.Orbit`
            Orbital elements.

        """

        desg = self.resolve_object(obj)[1]

        if isinstance(epochs, dict):
            _epochs = epochs
        else:
            _epochs = util.epochs_to_jd(epochs)

            d = np.diff(_epochs)
            if any(d <= 0):
                raise ValueError(
                    'Epoch dates must be increasing and unique: {}'.format(
                        _epochs))

            if len(_epochs) > 300:
                orb = None
                N = np.ceil(len(_epochs) / 200)
                for e in np.array_split(_epochs, N):
                    _orb = self.get_orbit_exact(obj, e, cache=cache)
                    if orb:
                        orb.add_rows(_orb)
                    else:
                        orb = _orb
                return orb

        kwargs = dict(epochs=_epochs, cache=cache)
        if Names.asteroid_or_comet(desg) == 'comet':
            kwargs.update(id_type='designation', closest_apparition=True,
                          no_fragments=True)

        orb = Orbit.from_horizons(desg, **kwargs)

        return orb

    def resolve_objects(self, objects):
        """Resolve objects to database object ID and designation.


        Parmeters
        ---------
        obj: array of str or int
            Objects to resolve: use strings for designation, ints for
            object ID.


        Returns
        -------
        objects : tuple
            (objid, desg) where object ID is ``None`` if not found.


        Raises
        ------
        ``BadObjectID`` if an object ID is provided but not in the
        database.

        """

        return tuple((self.resolve_object(obj) for obj in objects))

    def resolve_object(self, obj):
        """Resolve object to database object ID and designation.


        Parmeters
        ---------
        obj: str or int
            Object to resolve: use strings for designation, int for
            object ID.


        Returns
        -------
        objid: int
            Object ID or ``None`` if a designation is queried but not
            in the database.

        desg: str
            Object designation.


        Raises
        ------
        ``BadObjectID`` if an object ID is provided but not in the
        database.

        """
        if isinstance(obj, str):
            cmd = '''SELECT * FROM obj WHERE desg=?'''
            row = self.execute(cmd, [obj]).fetchone()
            if row is None:
                return None, str(obj)
            else:
                return int(row[0]), str(row[1])
        else:
            cmd = '''SELECT * FROM obj WHERE objid=?'''
            row = self.execute(cmd, [obj]).fetchone()
            if row is None:
                raise BadObjectID('{} not found in database'.format(obj))
            else:
                return int(row[0]), str(row[1])
