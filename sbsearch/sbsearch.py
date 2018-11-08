# Licensed with the 3-clause BSD license.  See LICENSE for details.
import numpy as np

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

    """

    def __init__(self, config=None, savelog=False):
        self.config = Config() if config is None else config
        fn = self.config['logfile'] if log else '/dev/null'
        self.logger = logging.setup(filename=fn)
        self.db = sqilte3.connect(config['filename'], factory=SBDB)
        self.db.verify_tables(self.config, self.logger)

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
            Integer step size in hours or days.  If `None`, then an
            adaptable step size will be used, based on distance to the
            observer.

        clean : bool, optional
            If ``True``, remove ephemerides, between and including
            start and stop, from the database.  With ``start is
            None``, remove all points before ``stop``.  With ``stop is
            None``, remove all points after ``start``.  If both are
            ``None``, all ephemeris points are removed.

        """

        objid = self.db.resolve_object(obj)[0]
        jd_start = None if start is None else Time(start, scale='utc').jd
        jd_stop = None if stop is None else Time(stop, scale='utc').jd
        if clean:
            n = self.db.clean_ephemeris(objid, jd_start, jd_stop)
            self.logger.info('{} rows deleted from eph table'.format(n))
        else:
            n = self.db.add_ephemeris(
                objid, self.config['location'], jd_start, jd_stop, step=step)
            self.logger.info('{} rows added to eph table'.format(n))

    def find_obs(self, obj, start=None, stop=None):
        """Find observations covering object.

        Parameters
        ----------
        objid : int or desg
            Object to find.

        start : float or `~astropy.time.Time`, optional
            Search after this date, inclusive.

        stop : float or `~astropy.time.Time`, optional
            Search before this date, inclusive.

        Returns
        -------
        obsids : array of int
            Observations with this object.

        """

        objid, desg = self.db.resolve_object(obj)
        columns = ('jd_start,jd_stop,ra,dec,ra1,dec1,ra2,dec2,'
                   'ra3,dec3,ra4,dec4')

        segments = self.db.get_ephemeris_segments(
            objid=objdid, start=start, stop=stop)

        found = []
        for segment in segments:
            obsids = self.db.get_observations_overlapping(box=segment)
            obs = list(zip(*self.db.get_observations_by_id(
                obsids, columns=columns)))

            epochs = (obs[0] + obs[1]) / 2
            eph = self.db.get_ephemeris_interp(objid, epochs)

            for i in range(len(obs)):
                corners = SkyCoord(obs[2::2], obs[3::2], 'rad')
                if interior_test(eph[i], corners):
                    found.append(obsids[i])

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
