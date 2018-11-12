# Licensed with the 3-clause BSD license.  See LICENSE for details.
import re
import sqlite3
import itertools
from logging import Logger

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
from .exceptions import BadObjectID


class SBDB(sqlite3.Connection):
    """Database object for SBSearch."""

    # observation table name
    OBS_TABLE = 'obs'

    # observation table columns; only append to this list
    OBS_COLUMNS = [
        'obsid INTEGER PRIMARY KEY',
        'jd_start FLOAT',
        'jd_stop FLOAT',
        'ra FLOAT',
        'dec FLOAT',
        'ra1 FLOAT',
        'dec1 FLOAT',
        'ra2 FLOAT',
        'dec2 FLOAT',
        'ra3 FLOAT',
        'dec3 FLOAT',
        'ra4 FLOAT',
        'dec4 FLOAT'
    ]

    def __init__(self, *args):
        super().__init__(*args)
        sqlite3.register_adapter(np.int64, int)
        sqlite3.register_adapter(np.int32, int)
        sqlite3.register_adapter(np.float64, float)
        sqlite3.register_adapter(np.float32, float)

        self.execute('PRAGMA foreign_keys = 1')
        self.execute('PRAGMA recursive_triggers = 1')
        self.row_factory = sqlite3.Row

    @classmethod
    def make_test_db(cls, target=':memory:'):
        db = sqlite3.connect(target, 5, 0, "DEFERRED", True, SBDB)
        db.verify_tables(Logger('test'))
        db.add_object('C/1995 O1')
        db.add_object('2P')

        # tile the sky
        N_tiles = 10
        sky_tiles = []
        ra_steps = np.linspace(0, 2 * np.pi, N_tiles + 1)
        dec_steps = np.linspace(-np.pi / 2, np.pi / 2, N_tiles + 1)
        sky_tiles = np.zeros((10, N_tiles**2))
        for i in range(N_tiles):
            for j in range(N_tiles):
                sky_tiles[:, i * N_tiles + j] = (
                    np.mean(ra_steps[i:i+2]),
                    np.mean(dec_steps[j:j+2]),
                    ra_steps[i], dec_steps[j],
                    ra_steps[i], dec_steps[j+1],
                    ra_steps[i+1], dec_steps[j+1],
                    ra_steps[i+1], dec_steps[j])
        obsids = range(N_tiles**2)
        start = 2458119.5 + np.arange(N_tiles**2) * 30 / 86400
        stop = start + 30 / 86400
        db.add_observations(columns=[obsids, start, stop] + list(sky_tiles))
        db.add_ephemeris(2, '500', 2458119.5, 2458121.5, step='1d',
                         source='jpl', cache=True)
        db.add_found(2, [1, 2, 3], '500', cache=True)
        return db

    def verify_tables(self, logger):
        """Verify SBSearch tables.

        Parameters
        ----------
        config : `~sbsearch.config.Config`
            Use this configuration set.

        logger : `~logging.Logger`
            Log messages to this logger.

        """

        c = self.execute("SELECT name FROM sqlite_master WHERE type='table'")
        names = list([row[0] for row in c])
        n = 0
        for name in ['obj', 'eph', 'eph_tree', self.OBS_TABLE,
                     self.OBS_TABLE + '_tree', self.OBS_TABLE + '_found']:
            if name in names:
                n += 1
            else:
                logger.error('table {} is missing'.format(name))

        if n != 6:
            self.create_tables(logger)

        logger.info('Tables verified.')

    def create_tables(self, logger):
        """Create sbsearch database tables.

        Inhereting classes will generally override this method to and
        call ``super().db_init`` with their own column definitions.

        Parameters
        ----------
        logger : `~logging.Logger`
            Log messages to this logger.

        Notes
        -----
        The following columns are required by ``SBSearch``:

            obsid INTEGER PRIMARY KEY
            jd_start FLOAT
            jd_stop FLOAT
            ra FLOAT
            dec FLOAT
            ra1 FLOAT
            dec1 FLOAT
            ra2 FLOAT
            dec2 FLOAT
            ra3 FLOAT
            dec3 FLOAT
            ra4 FLOAT
            dec4 FLOAT

        Where ra, dec is the center of the image, and raN, decN are
        the corners.

        """

        schema.create(self, self.OBS_TABLE, self.OBS_COLUMNS)
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
                groups = itertools.groupby(eph, lambda e: e['delta'] < limit)
                for inside, epochs in groups:
                    if not inside:
                        continue

                    jd = list([e['jd'] for e in epochs])
                    if len(jd) == 1:
                        # just one epoch inside delta limit, skip
                        continue

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
                    ra_rate = 'dRA cos(Dec)'
                    epochs = {'start': next_step, 'step': step, 'number': n}
                elif source == 'jpl':
                    ra_rate = 'dRA'
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
                          eph[ra_rate][i].to('arcsec/hr').value,
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

    def add_found(self, objid, obsids, location, update=False, cache=False):
        """Add found objects to found database.

        Parameters
        ----------
        objid : int
            Found object ID.

        obsids : array-like of int
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

        obsids = np.array(obsids)
        foundids = np.zeros(obsids.shape)

        # already in found database?
        if not update:
            missing = []
            for i in range(len(obsids)):
                f = self.execute('''
                SELECT foundid FROM {}_found
                WHERE objid=? AND obsid=?
                '''.format(self.OBS_TABLE), [objid, obsids[i]]).fetchone()
                if f:
                    foundids[i] = f[0]
                else:
                    missing.append(i)

            if len(missing) > 0:
                foundids[missing] = self.add_found(
                    objid, obsids[missing], location, update=True,
                    cache=cache)

            return foundids[missing]

        rows = self.get_observations_by_id(
            obsids, columns='jd_start,jd_stop', generator=True)
        jd = np.array(
            [(row['jd_start'] + row['jd_stop']) / 2 for row in rows])

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

        for i in range(len(obsids)):
            data = [objid, obsids[i], jd[i]]
            for col, unit in (('ra', 'deg'), ('dec', 'deg'),
                              ('dra', 'arcsec/hr'),
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
            INSERT OR REPLACE INTO {}_found
            VALUES (null,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
            '''.format(self.OBS_TABLE), data)
            foundids[i] = c.lastrowid

        self.commit()

        return foundids

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

    def add_observations(self, rows=None, columns=None, logger=None):
        """Add observations to observation table.

        RA, Dec must be in radians.  If observations already exist for
        a given observation ID, the new data are ignored.

        Parameters
        ----------
        rows : list or tuple of array-like
            Rows of data to insert.

        columns : list or tuple of array-like
            Columns of data to insert.

        logger : `~logging.Logger`, optional
            Report progress to this logger.


        """
        def row_iterator(rows, columns):
            if rows is not None:
                for row in rows:
                    yield row
            if columns is not None:
                for row in zip(*columns):
                    yield row

        obs_cmd = '''
        INSERT OR IGNORE INTO {} VALUES ({})
        '''.format(self.OBS_TABLE, ','.join('?' * len(self.OBS_COLUMNS)))

        tree_cmd = '''
        INSERT OR IGNORE INTO {}_tree VALUES (?,?,?,?,?,?,?,?,?)
        '''.format(self.OBS_TABLE)

        if logger is not None:
            tri = ProgressTriangle(1, logger, log=True)

        for row in row_iterator(rows, columns):
            c = self.execute(obs_cmd, row)
            obsid = c.lastrowid

            mjd = np.array((row[1], row[2])) - 2400000.5
            ra = [row[3], row[5], row[7], row[9], row[11]]
            dec = [row[4], row[6], row[8], row[10], row[12]]
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
            Object ID.

        jd_start, jd_stop: float or ``None``
            Julian date range(inclusive), or ``None`` for unbounded.

        Returns
        -------
        count: int
            Number of rows removed.

        """
        constraints = [('objid=?', objid)]
        constraints.extend(util.date_constraints(jd_start, jd_stop))
        cmd, parameters = util.assemble_sql(
            'DELETE FROM {}_found'.format(self.OBS_TABLE), [], constraints)

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
                            cache=False):
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
            Source to use: 'mpc' or 'jpl'.

        cache : bool, optional
            Use cached ephemerides; primarily for testing.

        Returns
        -------
        eph : `~sbpy.data.ephem.Ephem`
            Ephemeris.

        """

        if source not in ['mpc', 'jpl']:
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
            eph = Ephem.from_mpc(desg, _epochs, location=location,
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

        return eph

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

    def get_ephemeris_segments(self, objid=None, start=None, stop=None):
        """Get ephemeris segments.

        Parameters
        ----------
        objid : int, optional
            Find this object.

        start : float or `~astropy.time.Time`, optional
            Search starting at this epoch, inclusive.

        stop : float or `~astropy.time.Time`, optional
            Search upto and including at this epoch.

        Returns
        -------
        generator that returns ``(ephid, segment)``, where ``ephid``
        is an integer and ``segment`` is a dict.

        """

        cmd = 'SELECT eph_tree.ephid,mjd0,mjd1,x0,x1,y0,y1,z0,z1 FROM eph_tree'
        constraints = []
        parameters = []

        if objid is not None:
            cmd += ' INNER JOIN eph ON eph.ephid=eph_tree.ephid'
            constraints.append(('objid=?', objid))

        if start is not None:
            mjd = util.epochs_to_jd([start])[0] - 2400000.5
            constraints.append(('mjd1 >= ?', mjd))

        if stop is not None:
            mjd = util.epochs_to_jd([stop])[0] - 2400000.5
            constraints.append(('mjd0 <= ?', mjd))

        cmd, parameters = util.assemble_sql(cmd, parameters, constraints)

        c = self.execute(cmd, parameters)
        for row in util.iterate_over(c):
            segment = dict(row)
            del segment['ephid']
            yield row['ephid'], segment

    def get_found(self, obj=None, start=None, stop=None,
                  columns='*', generator=False):
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

        generator : bool, optional
            Return a generator rather a list.

        Returns
        -------
        rows : list or generator
            Found object table rows.

        """

        cmd = 'SELECT {} FROM {}_found'.format(columns, self.OBS_TABLE)
        constraints = util.date_constraints(start, stop, column='obsjd')
        if obj is not None:
            objid = self.resolve_object(obj)[0]
            constraints.append(('objid=?', objid))

        cmd, parameters = util.assemble_sql(cmd, [], constraints)
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

        cmd = 'SELECT {} FROM {}_found WHERE foundid=?'.format(
            columns, self.OBS_TABLE)

        def g(foundids):
            for foundid in foundids:
                r = self.execute(cmd, [foundid]).fetchone()
                if r is None:
                    raise StopIteration
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

        cmd = 'SELECT {} FROM {}_found WHERE obsid=?'.format(
            columns, self.OBS_TABLE)

        def g(obsids):
            for obsid in obsids:
                r = self.execute(cmd, [obsid]).fetchone()
                if r is None:
                    raise StopIteration
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

    def get_observations_by_id(self, obsids, columns='*', generator=False):
        """Get observations by observation ID.

        Parameters
        ----------
        obsids: array-like
            Observation IDs to retrieve.

        columns: string, optional
            Columns to retrieve, default all.

        generator: bool, optional
            Return a generator rather a list.

        Returns
        -------
        rows: list or generator
            Observation table rows.

        """

        cmd = 'SELECT {} FROM {} WHERE obsid=?'.format(
            columns, self.OBS_TABLE)

        def g(obsids):
            for obsid in obsids:
                r = self.execute(cmd, [obsid]).fetchone()
                if r is None:
                    raise StopIteration
                else:
                    yield r

        if generator:
            return g(obsids)
        else:
            return list(g(obsids))

    def get_observations_by_date(self, start, stop, columns='*',
                                 generator=False):
        """Get observations by observation date.

        Parameters
        ----------
        start, stop: string or `~astropy.time.Time`, optional
            Date range to search, inclusive.

        columns: string, optional
            Columns to retrieve, default all.

        generator: bool, optional
            Return a generator rather a list.

        Returns
        -------
        rows: list or generator
            Observation table rows.

        """

        obsids = self.get_observations_overlapping(
            start=start, stop=stop, generator=True)

        def g(obsids):
            for obsid in obsids:
                yield obsid[0]

        return self.get_observations_by_id(g(obsids), columns=columns,
                                           generator=generator)

    def get_observations_overlapping(self, box=None, ra=None, dec=None,
                                     start=None, stop=None, generator=False):
        """Find observations that overlap the given area / time.

        Parameters
        ----------
        box: dict-like, optional
            Parameters to test for overlap, keyed by column name.

        ra, dec: array-like float or `~astropy.coordinates.Angle`, optional
            Search this area.  Floats are radians.  Must be at least
            three points.

        start, stop: float or `~astropy.time.Time`, optional
            Search this time range.  Floats are Julian dates.

        generator: bool, optional
            Return a generator instead of a list.

        Returns
        -------
        obsid: list or generator
            Observation IDs that might overlap box.

        """

        if box is not None and any([c is not None for c in (ra, dec, start, stop)]):
            raise ValueError(
                'Only one of box or ra/dec/jd can be searched at once.')

        if ra is not None:
            if len(ra) < 3:
                raise ValueError('RA requires at least 3 points')

        if dec is not None:
            if len(dec) < 3:
                raise ValueError('Dec requires at least 3 points')

        cmd = 'SELECT obsid FROM {}_tree'.format(self.OBS_TABLE)
        constraints = []
        key2constraint = {
            'mjd0': 'mjd1 >= ?',
            'mjd1': 'mjd0 <= ?',
            'x0': 'x1 > ?',
            'x1': 'x0 < ?',
            'y0': 'y1 > ?',
            'y1': 'y0 < ?',
            'z0': 'z1 > ?',
            'z1': 'z0 < ?'
        }

        query = {}
        if box is not None:
            query.update(box)

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
            mjd = util.epochs_to_time([start]).mjd[0]
            query['mjd0'] = mjd
        if stop is not None:
            mjd = util.epochs_to_time([stop]).mjd[0]
            query['mjd1'] = mjd

        for k in query.keys():
            constraints.append((key2constraint[k], query[k]))

        cmd, parameters = util.assemble_sql(cmd, [], constraints)
        rows = self.execute(cmd, parameters)
        if generator:
            return util.iterate_over(rows)
        else:
            return list([row[0] for row in rows])

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
                eph = None
                N = np.ceil(len(_epochs) / 200)
                for e in np.array_split(_epochs, N):
                    _eph = self.get_orbit_exact(obj, e, cache=cache)
                    if eph:
                        eph.add_rows(_eph)
                    else:
                        eph = _eph
                return eph

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

    def test_observation_coverage(self, start, stop):
        """Test for any observations within the requested range.

        Parameters
        ----------
        start, stop: float or `~astropy.time.Time`, optional
            Search this time range.  Floats are Julian dates.
            ``None`` for unbounded limits.

        Returns
        -------
        coverage : bool

        """

        cmd = 'SELECT * FROM {}'.format(self.OBS_TABLE)
        constraints = util.date_constraints(start, stop, column='jd_start')
        cmd, parameters = util.assemble_sql(cmd, [], constraints)
        cmd += ' LIMIT 1'
        rows = self.execute(cmd, parameters).fetchone()
        return rows is not None
