# Licensed with the 3-clause BSD license.  See LICENSE for details.
import re
import sqlite3
import itertools

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
from sbpy.data import Ephem, Names

from . import util, schema


class SBDB(sqlite3.Connection):
    """Database object for SBSearch."""

    # observation table name
    obs_table = 'obs'

    # observation table columns; inhereting classes will generally
    # want to add to this list; the following are required:
    obs_columns = [
        'obsid INTEGER PRIMARY KEY',
        'obsjd FLOAT',
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
            obsjd FLOAT
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
                    objid, jd_start, jd_stop, columns=['jd', 'delta'])
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

                for i in range(len(eph)):
                    # save to ephemeris table
                    cursor = self.execute('''
                    INSERT OR REPLACE INTO eph
                    VALUES (null,?,?,?,?,?,?,?,?,?)
                    ''', (objid,
                          eph['Date'][i].jd,
                          eph['r'][i].value,
                          eph['Delta'][i].value,
                          eph['RA'][i].to('rad').value,
                          eph['Dec'][i].to('rad').value,
                          eph['dRA cos(Dec)'][i].value,
                          eph['ddec'][i].value,
                          today))

                    # save to ephemeris tree
                    indices = (max(0, i - 1), i, min(i, len(eph) - 1))
                    c = tuple((coords[j] for j in indices))
                    _jd = tuple((jd[j] for j in indices))
                    limits = util.eph_to_limits(c, _jd, half_step)
                    self.execute('''
                    INSERT OR REPLACE INTO eph_tree VALUES (
                      last_insert_rowid(),?,?,?,?,?,?,?,?
                    )''', limits)

                count += len(eph)
                next_step = jd_start + (step * (count + 1)).to(u.day).value

        return count

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

    def get_ephemeris(self, objid, jd_start, jd_stop, columns=None,
                      iterator=False, order=True):
        """Get ephemeris data from database.

        Parameters
        ----------
        objid : int
            Database object ID.

        jd_start, jd_stop : float
            Julian date range (inclusive).

        columns : list or tuple, optional
            List of columns to return.  Default all.

        iterator : bool, optional
            ``True`` to return an iterator, else return all rows.

        order : bool, optional
            ``True`` to sort by Julian date.

        Returns
        -------
        eph : iterator or list of rows

        """

        if columns is None:
            cmd = 'SELECT * FROM eph'
        else:
            if not isinstance(columns, (list, tuple)):
                raise ValueError('columns must be list or tuple')
            cmd = 'SELECT ' + ','.join(columns) + ' FROM eph'

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

    def get_ephemeris_exact(self, desg, location, epochs, source='mpc',
                            cache=False):
        """Generate ephemeris at specific epochs from external source.

        Parameters
        ----------
        desg : string
            Object designation.

        location : string
            Observer location.

        epochs : array-like or dict
            Compute ephemeris at these epochs.  For arrays, must be
            floats (for Julian date) or else parsable by
            `~astropy.time.Time`.  Dictionaries are passed to the
            ephemeris source.

        source : string, optional
            Source to use: 'mpc' or 'jpl'.

        cache : bool, optional
            Use cached ephemerides; primarily for testing.

        Returns
        -------
        eph : `~sbpy.data.ephem.Ephem`
            Ephemeris.

        Raises
        ------
        EphemerisError
            For bad ephemeris coverage at requested epochs.

        """

        if isinstance(epochs, dict):
            _epochs = epochs
        else:
            _epochs = util.epochs_to_time(epochs)

        if source == 'mpc':
            eph = Ephem.from_mpc(desg, epochs, location=location,
                                 proper_motion='sky',
                                 proper_motion_unit='rad/s',
                                 cache=cache)
        elif source == 'jpl':
            kwargs = dict(id_type='designation', epochs=_epochs,
                          quantities='1,3,9,19,20,23,24,36',
                          cache=cache)
            if Names.asteroid_or_comet(desg) == 'comet':
                kwargs.update(closest_apparition=True, no_fragments=True)

            eph = Ephem.from_horizons(desg, **kwargs)

        return eph

    def get_ephemeris_interp(self, obj, epochs):
        """Get ephemeris at specific epochs by interpolation.

        Parameters
        ----------
        obj : str or int
            Object designation or obsid.

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

        objid = self.object_id(obj)
        coords = []
        for epoch in util.epochs_to_time(epochs):
            # get all points within a day
            jd0 = epoch.jd
            rows = self._get_eph(objid, jd0 - 0.5, jd0 + 0.5)
            jd, ra, dec = list(zip(*rows))

            # find nearest 2 points
            i = np.abs(np.array(jd) - jd0).argmin()
            jd_i = jd.pop(i)
            c_i = SkyCoord(ra[i], dec[i], unit='deg')
            j = np.abs(np.array(jd) - jd0).argmin()
            jd_j = jd.pop(j)
            c_j = SkyCoord(ra[j], dec[j], unit='deg')
            del jd, ra, dec  # arrays no longer aligned

            c = util.spherical_interpolation(c_i, c_j, jd_i, jd_j, jd0)
            coords.append(c)

        return SkyCoord(coords)

    def clean_ephemeris(self, objid, jd_start, jd_stop):
        """Remove ephemeris between dates (inclusive).

        Parameters
        ----------
        obj : str or int
            Object designation or obsid.

        jd_start, jd_stop : float
            Julian date range (inclusive).

        Returns
        -------
        n : int
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
        obj : str or int
            Object to resolve: use strings for designation, int for
            object ID.

        Returns
        -------
        objid : int
            Object ID.

        desg : str
            Object designation.

        """
        if isinstance(obj, str):
            cmd = '''SELECT * FROM obj WHERE desg=?'''
        else:
            cmd = '''SELECT * FROM obj WHERE objid=?'''

        row = self.execute(cmd, [obj]).fetchone()
        return int(row[0]), str(row[1])
