# Licensed with the 3-clause BSD license.  See LICENSE for details.
import itertools
import sqlite3

import numpy as np
import astropy.units as u
from astroquery.mpc import MPC

from . import logging
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
        self._db_connect()

    def _db_connect(self):
        """Connect to sbsearch database.

        If none of the required databases exist, ``_db_init`` is
        called.  If any are missing, ``CorruptDB`` is raised.

        """
        sqlite3.register_adapter(np.int64, int)
        sqlite3.register_adapter(np.int32, int)
        sqlite3.register_adapter(np.float64, float)
        sqlite3.register_adapter(np.float32, float)

        filename = self.config['database']
        self.db = sqlite3.connect(filename)
        self.db.execute('PRAGMA foreign_keys = 1')
        self.db.execute('PRAGMA recursive_triggers = 1')
        self.db.row_factory = sqlite3.Row

        count = self.db.execute('''
        SELECT count() FROM sqlite_master
        WHERE type='table'
          AND (
            name='obj' OR name='eph' OR name='eph_tree' OR
            name=? OR name=? OR name=?
        )
        ''', [self.config.obs_table, self.config.obs_tree, self.config.found_table]
        ).fetchone()[0]
        if count != 6:
            self._db_init()

        self.logger.info('Connected to database: {}'.format(filename))

    def _db_init(self, obs_columns=None):
        """Create sbsearch database tables.

        Inhereting classes will generally override this method to and
        call ``super()._db_init`` with their own column definitions.


        Parameters
        ----------
        obs_columns : list of strings, optional
            Column definitions for the observation table.  See Notes
            for required columns.

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

        Where ra, dec are the center of the image, and raN, decN are
        the corners.

        """

        if obs_columns is None:
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

        self.db.executescript('''
CREATE TABLE IF NOT EXISTS obj(
  objid INTEGER PRIMARY KEY,
  desg TEXT
);

CREATE TABLE IF NOT EXISTS eph(
  ephid INTEGER PRIMARY KEY,
  objid INTEGER,
  jd FLOAT,
  rh FLOAT,
  delta FLOAT,
  ra FLOAT,
  dec FLOAT,
  dra FLOAT,
  ddec FLOAT,
  retrieved TEXT,
  FOREIGN KEY(objid) REFERENCES objects(objid)
);

CREATE VIRTUAL TABLE IF NOT EXISTS eph_tree USING RTREE(
  ephid INTEGER PRIMARY KEY,
  mjd0 FLOAT,
  mjd1 FLOAT,
  x0 FLOAT,
  x1 FLOAT,
  y0 FLOAT,
  y1 FLOAT,
  z0 FLOAT,
  z1 FLOAT
);

CREATE TRIGGER IF NOT EXISTS delete_eph_from_tree BEFORE DELETE ON eph
BEGIN
  DELETE FROM eph_tree WHERE ephid=old.ephid
END;

CREATE TRIGGER IF NOT EXISTS delete_object_from_eph BEFORE DELETE ON obj
BEGIN
  DELETE FROM eph WHERE objid=old.objid
END;

/* observation table */
CREATE TABLE IF NOT EXISTS ?(
  ''' + ',\n  '.join(obs_columns) + '''
);

/* observation rtree */
CREATE VIRTUAL TABLE IF NOT EXISTS ? USING RTREE(
  obsid INTEGER PRIMARY KEY,
  mjd0 FLOAT,
  mjd1 FLOAT,
  x0 FLOAT,
  x1 FLOAT,
  y0 FLOAT,
  y1 FLOAT,
  z0 FLOAT,
  z1 FLOAT
);

/* delete from obs tree before delete from obs table */
CREATE TRIGGER IF NOT EXISTS ? BEFORE DELETE ON ?
BEGIN
  DELETE FROM ? WHERE obsid=old.obsid
END;

/* objects found in obs table */
CREATE TABLE IF NOT EXISTS ?(
  foundid INTEGER PRIMARY KEY,
  objid INTEGER,
  obsid INTEGER,
  obsjd TEXT,
  ra FLOAT,
  dec FLOAT,
  dra FLOAT,
  ddec FLOAT,
  ra3sig FLOAT,
  dec3sig FLOAT,
  vmag FLOAT,
  rh FLOAT,
  rdot FLOAT,
  delta FLOAT,
  phase FLOAT,
  selong FLOAT,
  sangle FLOAT,
  vangle FLOAT,
  trueanomaly FLOAT,
  tmtp FLOAT,
  FOREIGN KEY(obsid) REFERENCES ?(obsid),
  FOREIGN KEY(objid) REFERENCES obj(objid)
);

/* delete from found before delete from obs */
CREATE TRIGGER IF NOT EXISTS ? BEFORE DELETE ON ?
BEGIN
  DELETE FROM ? WHERE obsid=old.obsid
END;

/* delete from found before delete from obj */
CREATE TRIGGER IF NOT EXISTS ? BEFORE DELETE ON ?
BEGIN
  DELETE FROM ? WHERE objid=old.objid
END;

''', (
            self.config.obs_table,
            self.config.obs_tree,
            'delete_obs_from_' + self.config.obs_tree,
            self.config.obs_table, self.config.obs_tree,
            self.config.found_table,
            self.config.obs_table,
            'delete_obs_from_' + self.config.found_table,
            self.config.obs_table, self.config.found_table,
            'delete_obj_from_' + self.config.found_table,
            self.config.obj_table, self.config.found_table
        ))
        self.logger.info('Created database tables and triggers.')

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
                [self.obs_table] + list(observation))
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
            Step size, in a format usable by
            `astroquery.mpc.MPC.get_ephemeris`.  If `None`, then an
            adaptable step size will be used, based on distance to the
            Earth.

        clean : bool, optional
            If ``True``, remove ephemerides, between and including
            start and stop, from the database.  With ``start is
            None``, remove all points before ``stop``.  With ``stop is
            None``, remove all points after ``start``.  If both are
            ``None``, all ephemeris points are removed.

        """

        objid = self.resolve_object(obj)
        jd_start = None if start is None else Time(start, scale='utc').jd
        jd_stop = None if stop is None else Time(stop, scale='utc').jd
        if clean:
            self._clean_eph(objid, jd_start, jd_stop)
        else:
            self._add_eph(objid, jd_start, jd_stop, step=step)

    def clean_found(self, objects, start=None, stop=None):
        """Remove found objects from the database.


        """

    def get_ephemeris(self, desg, epochs, exact=False):
        """Get ephemeris at ``epochs``."""
        pass

    def find_by_date(self, desg, dates, exact=True):
        """Find observations covering objects."""
        pass

    def find_in_obs(self, obsid, exact=True):
        """Find all objects in an observation."""
        pass

    def object_id(self, desg, add=False):
        """Translate designation into sbsearch object ID."""
        pass

    def interior_test(self, point, corners):
        """Test if point is interior to corners assuming spherical geometry."""

    def _add_eph(self, objid, jd_start, jd_stop, step=None):
        if jd_start is None:
            raise ValueError('Start epoch required for adding ephemerides.')
        if jd_stop is None:
            raise ValueError('Stop epoch required for adding ephemerides.')

        if step is None:
            # daily ephemeris for delta > 1
            self._add_eph(objid, jd_start, jd_stop, step='1d')

            # ZChecker analysis, Oct 2018: error ~ 2" / delta for 6h time step
            for limit, substep in zip((1, '4h'), (0.25, '1h')):
                eph = self._get_eph(
                    objid, jd_start, jd_stop, columns='jd,delta')
                groups = itertools.groupby(eph, lambda e: e['delta'] < limit)
                for inside, epochs in groups:
                    if not inside:
                        continue

                    jd = list([e['jd'] for e in epochs])
                    if len(jd) == 1:
                        continue

                    self._add_eph(objid, jd[0], jd[-1], step=substep)
        else:
            desg = self._get_desg(objid)
            step = u.Quantity(step)
            count = 0
            total = round((jd_stop - jd_start) / step.to('d').value + 1)
            next_step = jd_start
            mpc = MPC()
            today = Time.now().iso[:10]
            while count < total:
                n = total - count

                eph = mpc.get_ephemeris(
                    desg, location=self.config['location'], start=next_step,
                    step=step, number=n, proper_motion='sky', proper_motion_unit='rad/s',
                    cache=False)

                for i in range(len(eph)):
                    self.db.execute('''
                    INSERT OR REPLACE INTO eph VALUES (null,?,?,?,?,?,?,?,?,?)
                    ''', (objid, eph[i]['Date'].jd, eph[i]['r'], eph[i]['Delta'],
                          np.radians(eph[i]['RA']), np.radians(eph[i]['Dec']),
                          eph[i]['dRA cos(Dec)'], eph[i]['dDec'], today))

                count += len(eph)
                next_step = jd_start + (step * (count + 1)).to(u.day).value

    def _add_object(self, desg):
        pass

    def _get_desg(self, objid):
        pass

    def _get_eph(self, objid, jd_start, jd_stop, columns=None, iterator=False, order=True):
        if colums is None:
            cmd = 'SELECT * FROM eph WHERE'
        else:
            cmd = 'SELECT ' + columns + ' FROM eph WHERE'
        if objid is None:
            cmd += ' objid NOTNULL'
        else:
            cmd += ' objid=?'

        if jd_start is not None:
            cmd += ' AND jd>=?'
            parameters.append(jd_start)

        if jd_stop is not None:
            cmd += ' AND jd<=?'
            parameters.append(jd_stop)

        if cmd.endswith('WHERE'):
            cmd = cmd[:-5]

        if order:
            cmd += ' ORDER BY jd'

        c = self.db.execute(cmd, parameters)
        if iterator:
            self._iterate_over(c)
        else:
            return c.fetchall()

    def _clean_eph(self, objid, jd_start, jd_stop):
        cmd = 'DELETE FROM eph WHERE objid=?'
        parameters = [objid]

        if jd_start is not None:
            cmd += ' AND jd>=?'
            parameters.append(jd_start)

        if jd_stop is not None:
            cmd += ' AND jd<=?'
            parameters.append(jd_stop)

        c = self.db.execute(cmd, parameters)
        self.logger.info('{} rows deleted from eph table'.format(c.rowcount))

    def _iterate_over(self, cursor):
        while True:
            rows = cursor.fetchmany()
            if not rows:
                return
            for row in rows:
                yield row
