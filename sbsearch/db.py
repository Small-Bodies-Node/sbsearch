# Licensed with the 3-clause BSD license.  See LICENSE for details.
import re
import sqlite3
import itertools
from logging import Logger

import numpy as np
from numpy import pi
from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
from astropy.time import Time
from sbpy.data import Ephem, Names, Orbit

from . import util, schema


class SBDB(sqlite3.Connection):
    """Database object for SBSearch."""

    # observation table name
    obs_table = 'obs'

    # observation table columns; inhereting classes can append to this
    # list:
    obs_columns = [
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
    def make_test_db(cls):
        db = sqlite3.connect(':memory:', 5, 0, None, True, SBDB)
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
                         source='mpc', cache=True)
        db.add_found([1, 2, 3], [2, 2, 2], '500', cache=True)
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

        count = self.execute("""
        SELECT count() FROM sqlite_master
        WHERE type='table'
          AND (
            name='obj' OR name='eph' OR name='eph_tree' OR
            name='""" + self.obs_table + """' OR
            name='obs_tree' OR name='found_table'
          )""").fetchone()[0]
        if count != 6:
            self.create_tables(logger)
        else:
            logger.info('{} tables verified'.format(filename))

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

        schema.create(self, self.obs_table, self.obs_columns)
        logger.info('Created database tables and triggers.')

    def add_ephemeris(self, objid, location, jd_start, jd_stop, step=None,
                      source='mpc', cache=False):
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
                                       step='1d', source=source, cache=cache)

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
                        continue

                    count -= self.clean_ephemeris(objid, jd[0], jd[-1])
                    count += self.add_ephemeris(
                        objid, location, jd[0], jd[-1], step=substep,
                        source=source, cache=cache)
        else:
            desg = self.resolve_object(objid)[1]
            step = u.Quantity(step)
            count = 0
            total = round((jd_stop - jd_start) / step.to('day').value + 1)
            next_step = jd_start
            today = Time.now().iso[:10]
            while count < total:
                n = total - count

                epochs = {'start': next_step, 'step': step, 'number': n}
                eph = self.get_ephemeris_exact(
                    desg, location, epochs, source=source, cache=cache)

                jd = eph['Date'].jd
                coords = SkyCoord(eph['RA'], eph['Dec'])
                half_step = step.to('d').value / 2

                vmag = util.vmag_from_eph(eph)

                for i in range(len(eph)):
                    # save to ephemeris table
                    cursor = self.execute('''
                    INSERT OR REPLACE INTO eph
                    VALUES (null,?,?,?,?,?,?,?,?,?,?)
                    ''', (objid,
                          eph['Date'][i].jd,
                          eph['r'][i].value,
                          eph['Delta'][i].value,
                          eph['RA'][i].to('rad').value,
                          eph['Dec'][i].to('rad').value,
                          eph['dRA cos(Dec)'][i].value,
                          eph['ddec'][i].value,
                          vmag[i],
                          today))

                    # save to ephemeris tree
                    indices = (max(0, i - 1), i, min(i, len(eph) - 1))
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

    def add_found(self, obsids, objids, location, cache=False):
        """Add found objects to found database.

        Parameters
        ----------
        obsids : array-like of int
            Observations with found objects.

        objids : array-like of int
            Found object IDs.

        location : string
            Observer location.

        cache : bool, optional
            Use cached ephemerides; primarily for testing.

        Returns
        -------
        foundids : list
            Generated found IDs.

        """

        obsids = np.array(obsids)
        objids = np.array(objids)
        foundids = np.empty(obsids.shape)

        rows = self.get_observations_by_id(
            obsids, columns='jd_start,jd_stop', generator=True)
        jd = np.array(
            [(row['jd_start'] + row['jd_stop']) / 2 for row in rows])
        epochs = Time(jd, format='jd')

        for objid in np.unique(objids):
            i = np.flatnonzero(objids == objid)
            eph = self.get_ephemeris_exact(objid, location, epochs[i],
                                           source='jpl', cache=cache)
            orb = self.get_orbit_exact(objid, epochs[i], cache=cache)

            vmag = util.vmag_from_eph(eph)
            sangle = Angle(eph['sunTargetPA'] - 180 * u.deg)
            sangle = sangle.wrap_at(360 * u.deg).deg
            vangle = Angle(eph['velocityPA'] - 180 * u.deg)
            vangle = vangle.wrap_at(360 * u.deg).deg
            Tp = Time(orb['Tp_jd'], format='jd', scale='tt').utc.jd
            tmtp = Tp - jd[i]

            for j, k in enumerate(i):
                data = [objid, obsids[k], jd[k]]
                for col, unit in (('ra', 'deg'), ('dec', 'deg'),
                                  ('dra', 'arcsec/hr'),
                                  ('ddec', 'arcsec/hr'),
                                  ('RA_3sigma', 'arcsec'),
                                  ('DEC_3sigma', 'arcsec')):
                    data.append(eph[col][j].to(unit).value)
                data.append(vmag[j])
                for col, unit in (('r', 'au'), ('r_rate', 'km/s'),
                                  ('delta', 'au'), ('alpha', 'deg'),
                                  ('elong', 'deg')):
                    data.append(eph[col][j].value)
                data.extend((sangle[j], vangle[j], orb['nu'][j].value,
                             tmtp[j]))

                c = self.execute('''
                INSERT OR REPLACE INTO {}_found
                VALUES (null,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
                '''.format(self.obs_table), data)
                foundids[k] = c.lastrowid

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

    def add_observations(self, rows=None, columns=None):
        """Add observations to observation table.

        Parameters
        ----------
        rows : list or tuple of array-like
            Rows of data to insert.

        columns : list or tuple of array-like
            Columns of data to insert.

        """
        def row_iterator(rows, columns):
            if rows is not None:
                for row in rows:
                    yield row
            if columns is not None:
                for row in zip(*columns):
                    yield row

        obs_cmd = '''
        INSERT OR REPLACE INTO {} VALUES ({})
        '''.format(self.obs_table, ','.join('?' * len(self.obs_columns)))

        tree_cmd = '''
        INSERT OR REPLACE INTO {}_tree VALUES (?,?,?,?,?,?,?,?,?)
        '''.format(self.obs_table)

        for row in row_iterator(rows, columns):
            c = self.execute(obs_cmd, row)
            obsid = c.lastrowid

            mjd = np.array((row[1], row[2])) - 2450000.5
            ra = [row[3], row[5], row[7], row[9], row[11]]
            dec = [row[4], row[6], row[8], row[10], row[12]]
            xyz = util.rd2xyz(ra, dec)
            self.execute(
                tree_cmd, [obsid, min(mjd), max(mjd),
                           min(xyz[0]), max(xyz[0]),
                           min(xyz[1]), max(xyz[1]),
                           min(xyz[2]), max(xyz[2])])

    def get_ephemeris(self, objid, jd_start, jd_stop, columns='*',
                      iterator=False, order=True):
        """Get ephemeris data from database.

        Parameters
        ----------
        objid : int
            Database object ID.

        jd_start, jd_stop : float
            Julian date range (inclusive).

        columns : string, optional
            Columns to return.  Default all.

        iterator : bool, optional
            ``True`` to return an iterator, else return all rows.

        order : bool, optional
            ``True`` to sort by Julian date.

        Returns
        -------
        eph : iterator or list of rows

        """

        cmd = 'SELECT {} FROM eph'.format(columns)

        constraints = []
        parameters = []
        if objid is None:
            constraints.append(('objid ?', 'NOTNULL'))
        else:
            constraints.append(('objid=?', objid))

        constraints.extend(util.date_constraints(jd_start, jd_stop))
        cmd, parameters = util.assemble_sql(cmd, parameters, constraints)

        if order:
            cmd += ' ORDER BY jd'

        c = self.execute(cmd, parameters)
        if iterator:
            util.iterate_over(c)
        else:
            return c.fetchall()

    def get_ephemeris_exact(self, obj, location, epochs, source='mpc',
                            cache=False):
        """Generate ephemeris at specific epochs from external source.

        Parameters
        ----------
        obj : string or int
            Object designation or object ID.

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

        desg = self.resolve_object(obj)[1]

        if isinstance(epochs, dict):
            _epochs = epochs
        elif isinstance(epochs, Time):
            _epochs = epochs.jd
        else:
            _epochs = util.epochs_to_time(epochs).jd

        if source == 'mpc':
            eph = Ephem.from_mpc(desg, epochs, location=location,
                                 proper_motion='sky',
                                 proper_motion_unit='rad/s',
                                 cache=cache)
        elif source == 'jpl':
            kwargs = dict(id_type='designation', epochs=_epochs,
                          quantities='1,3,8,9,19,20,23,24,27,36',
                          cache=cache)
            if Names.asteroid_or_comet(desg) == 'comet':
                kwargs.update(closest_apparition=True, no_fragments=True)

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
        coords : `~astropy.coordinates.SkyCoord`
            RA and Dec at each epoch.

        Raises
        ------
        EphemerisError
            For bad ephemeris coverage at requested epochs.

        """

        coords = []
        for epoch in util.epochs_to_time(epochs):
            # get two nearest points to epoch
            jd0 = epoch.jd
            rows = self.execute('''
            SELECT jd, ra, dec, ABS(jd - ?) AS dt FROM eph
            WHERE objid=?
              AND dt < 5
            ORDER BY dt LIMIT 2
            ''', [jd0, objid])

            jd, ra, dec, dt = list(zip(*rows))

            i = np.argmin(jd)
            j = np.argmax(jd)

            a = SkyCoord(ra[i], dec[i], unit='deg')
            b = SkyCoord(ra[j], dec[j], unit='deg')

            c = util.spherical_interpolation(a, b, jd[i], jd[j], jd0)
            coords.append(c)

        return SkyCoord(coords)

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
        rows : iterable
            Rows from eph_tree table.

        """

        cmd = 'SELECT eph_tree.ephid,mjd0,mjd1,x0,x1,y0,y1,z0,z1 FROM eph_tree'
        constraints = []
        parameters = []

        if objid is not None:
            cmd += ' INNER JOIN eph ON eph.ephid=eph_tree.ephid'
            constraints.append(('objid=?', objid))

        if start is not None:
            mjd = util.epochs_to_time([start]).jd[0] - 2450000.5
            constraints.append(('mjd1 >= ?', mjd))

        if stop is not None:
            mjd = util.epochs_to_time([stop]).jd[0] - 2450000.5
            constraints.append(('mjd0 <= ?', mjd))

        cmd, parameters = util.assemble_sql(cmd, parameters, constraints)

        c = self.execute(cmd, parameters)
        return util.iterate_over(c)

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
            columns, self.obs_table)

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

    def get_found_by_date(self, start, stop, columns='*', generator=False):
        """Get found objects by observation dates.

        Parameters
        ----------
        start, stop : string or `~astropy.time.Time`, optional
            Date range to search, inclusive.

        columns : string, optional
            Columns to retrieve, default all.

        generator : bool, optional
            Return a generator rather a list.

        Returns
        -------
        rows : list or generator
            Found object table rows.

        """

        cmd = 'SELECT {} FROM {}_found WHERE obsjd>=? AND obsjd<=?'.format(
            columns, self.obs_table)
        rows = self.execute(cmd, util.epochs_to_time([start, stop]).jd)

        if generator:
            return util.iterate_over(rows)
        else:
            return list(rows.fetchall())

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
            columns, self.obs_table)

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

    def get_found_by_object(self, obj, columns='*', generator=False):
        """Get found objects by object name/ID.

        Parameters
        ----------
        obj : int or string
            Object to find.

        columns : string, optional
            Columns to retrieve, default all.

        generator : bool, optional
            Return a generator rather a list.

        Returns
        -------
        rows : list or generator
            Found object table rows.

        """

        cmd = 'SELECT {} FROM {}_found WHERE objid=?'.format(
            columns, self.obs_table)
        rows = self.execute(cmd, [self.resolve_object(obj)[0]])

        if generator:
            return util.iterate_over(rows)
        else:
            return list(rows.fetchall())

    def get_observations_by_id(self, obsids, columns='*', generator=False):
        """Get observations by observation ID.

        Parameters
        ----------
        obsids : array-like
            Observation IDs to retrieve.

        columns : string, optional
            Columns to retrieve, default all.

        generator : bool, optional
            Return a generator rather a list.

        Returns
        -------
        rows : list or generator
            Observation table rows.

        """

        cmd = 'SELECT {} FROM {} WHERE obsid=?'.format(
            columns, self.obs_table)

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
        start, stop : string or `~astropy.time.Time`, optional
            Date range to search, inclusive.

        columns : string, optional
            Columns to retrieve, default all.

        generator : bool, optional
            Return a generator rather a list.

        Returns
        -------
        rows : list or generator
            Observation table rows.

        """

        obsids = self.get_observations_overlapping(
            epochs=[start, stop], generator=True)

        def g(obsids):
            for obsid in obsids:
                yield obsid[0]

        return self.get_observations_by_id(g(obsids), columns=columns,
                                           generator=generator)

    def get_observations_overlapping(self, box=None, ra=None, dec=None,
                                     epochs=None, generator=False):
        """Find observations that overlap the given area / time.

        Parameters
        ----------
        box: dict-like, optional
            Parameters to test for overlap, keyed by column name.

        ra, dec: array-like float or `~astropy.coordinates.Angle`, optional
            Search this area.  Floats are radians.  Must be at least
            three points.

        epochs: iterable of float or `~astropy.time.Time`, optional
            Search this time range.  Floats are Julian dates.

        generator : bool, optional
            Return a generator instead of a list.

        Returns
        -------
        obsid: list or generator
            Observation IDs that might overlap box.

        """

        if box is not None and any([c is not None for c in (ra, dec, jd)]):
            raise ValueError(
                'Only one of box or ra/dec/jd can be searched at once.')

        cmd = 'SELECT obsid FROM {}_tree'.format(self.obs_table)
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

        box = {}
        if ra is not None and dec is None:
            if len(ra) < 3:
                raise ValueError('RA requires at least 3 points')
            cra = np.cos(ra)
            sra = np.sin(ra)
            box['x0'] = min(cra)
            box['x1'] = max(cra)
            box['y0'] = min(sra)
            box['y1'] = max(sra)
            box['z0'] = -1
            box['z1'] = 1
        elif ra is None and dec is not None:
            if len(dec) < 3:
                raise ValueError('Dec requires at least 3 points')
            sdec = np.sin(dec)
            box['x0'] = -1
            box['x1'] = 1
            box['y0'] = -1
            box['y1'] = 1
            box['z0'] = min(sdec)
            box['z1'] = max(sdec)
        elif ra is not None and dec is not None:
            if len(ra) < 3:
                raise ValueError('RA requires at least 3 points')
            if len(dec) < 3:
                raise ValueError('Dec requires at least 3 points')
            if len(ra) != len(dec):
                raise ValueError('RA and Dec must have the same length.')
            x, y, z = util.rd2xyz(ra, dec)
            box['x0'] = min(x)
            box['x1'] = max(x)
            box['y0'] = min(y)
            box['y1'] = max(y)
            box['z0'] = min(z)
            box['z1'] = max(z)
        if epochs is not None:
            _jd = util.epochs_to_time(epochs).jd
            box['mjd0'] = min(_jd) - 2450000.5
            box['mjd1'] = max(_jd) - 2450000.5

        for k in box.keys():
            constraints.append((key2constraint[k], box[k]))

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
        obj : string or int
            Object designation or object ID.

        epochs : array-like or dict
            Compute ephemeris at these epochs.  For arrays, must be
            floats (for Julian date) or else parsable by
            `~astropy.time.Time`.  Dictionaries are passed to the
            ephemeris source as is.

        cache : bool, optional
            Use cached ephemerides; primarily for testing.

        Returns
        -------
        orb : `~sbpy.data.Orbit`
            Orbital elements.

        """

        desg = self.resolve_object(obj)[1]

        if isinstance(epochs, dict):
            _epochs = epochs
        elif isinstance(epochs, Time):
            _epochs = epochs.jd
        else:
            _epochs = util.epochs_to_time(epochs).jd

        kwargs = dict(id_type='designation', epochs=_epochs,
                      cache=cache)
        if Names.asteroid_or_comet(desg) == 'comet':
            kwargs.update(closest_apparition=True, no_fragments=True)

        orb = Orbit.from_horizons(desg, **kwargs)

        return orb

    def clean_ephemeris(self, objid, jd_start, jd_stop):
        """Remove ephemeris between dates(inclusive).

        Parameters
        ----------
        obj: str or int
            Object designation or obsid.

        jd_start, jd_stop: float
            Julian date range(inclusive).

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

    def resolve_object(self, obj):
        """Resolved object to database object ID and designation.

        Parmeters
        ---------
        obj: str or int
            Object to resolve: use strings for designation, int for
            object ID.

        Returns
        -------
        objid: int
            Object ID.

        desg: str
            Object designation.

        """
        if isinstance(obj, str):
            cmd = '''SELECT * FROM obj WHERE desg=?'''
        else:
            cmd = '''SELECT * FROM obj WHERE objid=?'''

        row = self.execute(cmd, [obj]).fetchone()
        return int(row[0]), str(row[1])
