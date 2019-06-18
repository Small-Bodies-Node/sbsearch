# Licensed with the 3-clause BSD license.  See LICENSE for details.
from logging import Logger

import numpy as np

import sqlalchemy as sa
from sqlalchemy.orm.exc import NoResultFound

from astropy.coordinates import Angle
import astropy.units as u
from astropy.time import Time
from astropy.table import vstack
from sbpy.data import Ephem

from . import util, schema, ephem
from .schema import Obj, Eph, Found, Obs, GenericObs
from .util import RADec
from .logging import ProgressTriangle
from .exceptions import (
    BadObjectID,
    NoEphemerisError,
    SourceNotFoundError,
)


class SBDB:
    """Database object for SBSearch."""

    DB_NAMES = ['obj', 'eph', 'obs', 'generic_obs', 'found']

    def __init__(self, url, *args):
        self.engine = sa.create_engine(url, *args)
        Session = sa.orm.sessionmaker(bind=self.engine)
        self.session = Session()

    def __del__(self):
        self.close()

    def close(self):
        self.session.close()

    @classmethod
    def create_test_db(cls):
        """Create a test database.

        Requires database named "sbsearch_test".

        Scheme for handling temporary databases from:
        http://alextechrants.blogspot.com/2014/01/unit-testing-sqlalchemy-apps-part-2.html

        """
        from .test.skytiles import sky_tiles, N_tiles

        url = "postgresql:///sbsearch_test"
        db = SBDB(url)

        # remove previous temp work
        metadata = sa.MetaData(bind=db.engine)
        metadata.reflect()
        for table in metadata.tables.keys():
            if table == 'spatial_ref_sys':
                continue
            db.engine.execute('DROP TABLE IF EXISTS {} CASCADE'
                              .format(table))

        # create test database
        db.verify_database(Logger('test'))
        db.add_object('C/1995 O1')
        db.add_object('2P')

        # 30 s exposures
        exptime = 30 / 86400
        for i in range(N_tiles**2):
            obs = GenericObs(
                obsid=i,
                jd_start=2458119.5 + exptime * i,
                jd_stop=2458119.5 + exptime * (i + 1),
                fov=sky_tiles[i],
                filter='r',
                seeing=1.5,
                airmass=1.3,
                maglimit=25)
            db.session.add(obs)

        db.add_ephemeris(2, '500', 2458119.5, 2458121.5, step='1d',
                         source='jpl', cache=True)
        db.add_found_by_id(2, [1, 2, 3], '500', cache=True)
        return db

    def verify_database(self, logger, names=[]):
        """Verify SBSearch tables, triggers, etc.

        Parameters
        ----------
        logger : `~logging.Logger`
            Log messages to this logger.

        names : list, optional
            Additional database names to consider.

        """

        conn = self.engine.connect()
        metadata = sa.MetaData()
        metadata.reflect(self.engine)

        missing = False
        for name in self.DB_NAMES + names:
            if name not in metadata.tables.keys():
                missing = True
                logger.error('{} is missing from database'.format(name))

        if missing:
            schema.create(self.engine)
            logger.info('Created database tables and triggers.')

    def add_ephemeris(self, objid, location, start, stop, step=None,
                      source='jpl', cache=False):
        """Add ephemeris data to databse.

        Parameters
        ----------
        objid : int
            Database object ID.

        location : string
            Observer location.

        start, stop : float, string, `~astropy.time.Time`
            Observation date range (inclusive).

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

        desg = self.resolve_object(objid)[1]
        today = Time.now().iso[:10]
        count = 0

        t = util.epochs_to_time((start, stop))
        epochs = {
            'start': t[0].iso,
            'stop': t[1].iso,
            'step': step
        }
        eph = ephem.generate(desg, location, epochs, source=source,
                             cache=cache)

        jd = util.epochs_to_jd(eph['Date'].value)
        coords = RADec(eph['RA'], eph['Dec'])
        vmag = util.vmag_from_eph(eph)

        for i in range(len(eph)):
            # save ephemeris segment, which is used for
            # searching, yes they overlap
            p0 = max(i - 1, 0)
            p1 = min(len(eph) - 1, i + 1)

            segment = 'SRID=40001;LINESTRING({} {}, {} {})'.format(
                eph['RA'][p0].to('deg').value,
                eph['Dec'][p0].to('deg').value,
                eph['RA'][p1].to('deg').value,
                eph['Dec'][p1].to('deg').value)

            # ephemeris table data
            data = Eph(
                objid=objid,
                jd=jd[i],
                rh=eph['r'][i].value,
                delta=eph['Delta'][i].value,
                ra=eph['RA'][i].to('rad').value,
                dec=eph['Dec'][i].to('rad').value,
                dra=eph['dRA cos(Dec)'][i].to('arcsec/hr').value,
                ddec=eph['ddec'][i].to('arcsec/hr').value,
                unc_a=eph['SMAA_3sigma'][i].to('rad').value,
                unc_b=eph['SMIA_3sigma'][i].to('rad').value,
                unc_theta=eph['Theta_3sigma'][i].to('rad').value,
                vmag=vmag[i],
                segment=segment,
                retrieved=today
            )
            self.session.add(data)
            count += 1

        self.session.commit()

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

        # already in found database?
        if not update:
            foundids = np.zeros(len(observations), int)
            missing = []
            obs_missing = []
            for i, obs in enumerate(observations):
                try:
                    foundids[i] = (self.session.query(Found.foundid)
                                   .filter(Found.objid == objid)
                                   .filter(Found.obsid == obs.obsid)
                                   .one())[0]
                except NoResultFound:
                    missing.append(i)
                    obs_missing.append(obs)

            if len(missing) > 0:
                foundids[missing] = self.add_found(
                    objid, obs_missing, location, update=True,
                    cache=cache)

            return foundids[missing]

        jd = np.array([(obs.jd_start + obs.jd_stop) / 2
                       for obs in observations])

        jd_sorted, unsort_jd = np.unique(jd, return_inverse=True)
        desg = self.resolve_object(objid)[1]
        eph = ephem.generate(desg, location, jd_sorted, source='jpl',
                             cache=cache)
        orb = ephem.generate_orbit(desg, jd_sorted, cache=cache)

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

        found = []
        for i, obs in enumerate(observations):
            found.append(Found(
                objid=objid,
                obsid=obs.obsid,
                jd=jd[i],
                ra=eph['ra'][i].to('deg').value,
                dec=eph['dec'][i].to('deg').value,
                dra=eph['dRA cos(Dec)'][i].to('arcsec/hr').value,
                ddec=eph['ddec'][i].to('arcsec/hr').value,
                unc_a=eph['SMAA_3sigma'][i].to('arcsec').value,
                unc_b=eph['SMIA_3sigma'][i].to('arcsec').value,
                unc_theta=eph['Theta_3sigma'][i].to('deg').value,
                vmag=vmag[i],
                rh=eph['r'][i].to('au').value,
                rdot=eph['r_rate'][i].to('km/s').value,
                delta=eph['Delta'][i].to('au').value,
                phase=eph['alpha'][i].to('deg').value,
                selong=eph['elong'][i].to('deg').value,
                sangle=sangle[i],
                vangle=vangle[i],
                trueanomaly=orb['nu'][i].to('deg').value,
                tmtp=tmtp[i]
            ))

            self.session.add(found[-1])

        self.session.commit()

        foundids = [f.foundid for f in found]
        return foundids

    def add_found_by_id(self, objid, obsids, location, **kwargs):
        """Add found objects to found database using observation ID.

        Parameters
        ----------
        objid : int
            Found object ID.

        obsids : array-like
            Observation IDs with found object.

        location : string
            Observer location.

        **kwargs
            Any ``add_found`` keyword arguments.

        Returns
        -------
        foundids : list
            New found IDs.  If ``update`` is ``False``, found IDs that
            already exist will not be returned.

        """

        observations = self.get_observations_by_id(obsids)
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

        obj = Obj(desg=desg)
        self.session.add(obj)
        self.session.commit()
        return obj.objid

    def add_observations(self, observations, logger=None):
        """Add observations to database.

        RA, Dec must be in radians.  If observations already exist for
        a given observation ID, the new data are ignored.


        Parameters
        ----------
        observations : list of Obs
            Observations to insert.

        logger : `~logging.Logger`, optional
            Report progress to this logger.

        """

        if logger is not None:
            tri = ProgressTriangle(1, logger, base=10)

        # autoincrement work around for Postgres
        if self.engine.dialect.name == 'postgresql':
            self.session.execute('''
            SELECT setval('obs_obsid_seq', MAX(obsid)) FROM obs 
            ''')

        for obs in observations:
            self.session.add(obs)

            if logger is not None:
                tri.update()

        self.session.commit()

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

        eph = self.session.query(Eph).filter_by(objid=objid)
        eph = util.filter_by_date_range(eph, jd_start, jd_stop, Eph.jd).all()

        count = len(eph)
        self.session.delete(eph)
        self.session.commit()

        return count

    def clean_found(self, objid, jd_start, jd_stop):
        """Remove found objects.

        Parameters
        ----------
        objid : int
            Object ID or ``None`` for all objects.

        jd_start, jd_stop : float or ``None``
            Julian date range(inclusive), or ``None`` for unbounded.

        Returns
        -------
        count : int
            Number of rows removed.

        """

        found = self.session.query(Found).filter_by(objid=objid)
        found = util.filter_by_date_range(
            found, jd_start, jd_stop, Found.jd)

        count = 0
        for f in found:
            count += 1
            self.session.delete(f)
        self.session.commit()

        return count

    def get_ephemeris(self, objid, jd_start, jd_stop):
        """Get ephemeris data from database.

        Parameters
        ----------
        objid : int or ``None``
            Database object ID or ``None`` for any object.

        jd_start, jd_stop : float
            Julian date range (inclusive).

        Returns
        -------
        eph : sqlalchemy Query

        """

        eph = self.session.query(Eph)

        if objid is not None:
            eph = eph.filter_by(objid=objid)

        eph = util.filter_by_date_range(eph, jd_start, jd_stop, Eph.jd)

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

        eph = self.session.query(sa.func.min(Eph.jd).label('jd_min'),
                                 sa.func.max(Eph.jd).label('jd_max'))

        if objids:
            eph = eph.filter(Eph.objid.in_(objids))

        jd_min, jd_max = eph.one()
        if None in [jd_min, jd_max]:
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
            dt = sa.func.abs(Eph.jd - jd0).label('dt')
            eph = (
                self.session.query(Eph, dt)
                .filter(Eph.objid == objid)
                .filter(dt < 5)
                .order_by(dt)
                .limit(2)
                .all()
            )

            jd = [row[0].jd for row in eph]
            _ra = [row[0].ra for row in eph]
            _dec = [row[0].dec for row in eph]
            _vmag = [row[0].vmag for row in eph]
            dt = [row[1] for row in eph]

            i = np.argmin(jd)
            j = np.argmax(jd)

            a = RADec(_ra[i], _dec[i], unit='rad')
            b = RADec(_ra[j], _dec[j], unit='rad')

            c = util.spherical_interpolation(a, b, jd[i], jd[j], jd0)
            ra.append(c.ra)
            dec.append(c.dec)

            vmag.append(np.interp(jd0, jd, _vmag))

        return RADec(ra, dec, unit='rad'), np.array(vmag)

    def get_found(self, obj=None, start=None, stop=None):
        """Get found objects by object and/or date range.

        Parameters
        ----------
        obj : int or string, optional
            Find detections of this object.

        start, stop : string or `~astropy.time.Time`, optional
            Date range to search, inclusive.  ``None`` for unbounded
            limit.

        Returns
        -------
        found : sqalchemy Query
            Matching found objects.

        """

        start, stop = util.epochs_to_jd((start, stop))

        found = self.session.query(Found)

        if obj is not None:
            objid = self.resolve_object(obj)[0]
            found = found.join(Obj).filter(Found.objid == objid)

        if start is not None:
            found = found.filter(Found.jd >= start)

        if stop is not None:
            found = found.filter(Found.jd <= stop)

        return found

    def get_found_by_id(self, foundids):
        """Get found objects by found ID.

        Parameters
        ----------
        foundids : array-like, optional
            Found IDs to retrieve.

        Returns
        -------
        found : sqalchemy Query
            Matching found objects.

        """

        found = self.session.query(Found).filter(
            Found.foundid.in_(foundids))
        return found

    def get_found_by_obsid(self, obsids):
        """Get found objects by observation ID.

        Parameters
        ----------
        obsids : array-like
            Observation IDs to search.

        Returns
        -------
        found : sqalchemy Query
            Matching found objects.

        """

        found = self.session.query(Found).filter(
            Found.objid.in_(obsids))
        return found

    def get_objects(self):
        """Return all objects.

        Returns
        -------
        objects : sqlalchemy Query
            All objects.

        """

        return self.session.query(Obj)

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

        query = self.session.query(
            sa.func.min(Obs.jd_start),
            sa.func.max(Obs.jd_stop))

        if source:
            query = query.filter(Obs.source == source)

        jd_min, jd_max = query.one()
        if None in [jd_min, jd_max]:
            if source:
                msg = 'No observations for source: ' + source
            else:
                msg = 'No observations in database'
            raise SourceNotFoundError(msg)

        return jd_min, jd_max

    def get_observations_by_id(self, obsids):
        """Get observations by observation ID.

        Parameters
        ----------
        obsids: array-like
            Observation IDs to retrieve.

        Returns
        -------
        obs : sqlalchemy Query
            Matched observations.

        """

        obs = self.session.query(Obs).filter(Obs.obsid.in_(obsids))
        return obs

    def get_observations_by_date(self, start=None, stop=None):
        """Get observations by observation date.

        Parameters
        ----------
        start, stop: string or `~astropy.time.Time`, optional
            Date range to search, inclusive.

        Returns
        -------
        obs : sqlalchemy Query
            Matched observations.

        """

        start, stop = util.epochs_to_jd((start, stop))

        obs = self.session.query(Obs)

        if start is not None:
            obs = obs.filter(Obs.jd_stop >= start)

        if stop is not None:
            obs = obs.filter(Obs.jd_start <= stop)

        return obs

    def get_observations_covering(self, shape, start=None, stop=None):
        """Find observations covering the given shape.

        Parameters
        ----------
        shape : string
            PostGIS WKT spatial object.

        start, stop: float or `~astropy.time.Time`, optional
            Search this time range.  Floats are Julian dates.

        Returns
        -------
        obs : sqlalchemy Query
            Matching observations.

        """

        if start is not None or stop is not None:
            obs = self.get_observations_by_date(start, stop)

        obs = obs.filter(Obs.fov.ST_Covers(shape))

        return obs

    def get_observations_intersecting(self, shape, start=None, stop=None):
        """Find observations intersecting the given shape.

        Parameters
        ----------
        shape : string
            PostGIS WKT spatial object.

        start, stop: float or `~astropy.time.Time`, optional
            Search this time range.  Floats are Julian dates.

        Returns
        -------
        obs : sqlalchemy Query
            Matching observations.

        """

        if start is not None or stop is not None:
            obs = self.get_observations_by_date(start, stop)

        obs = obs.filter(Obs.fov.ST_Intersects(shape))

        return obs

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

        query = self.session.query(Obj.objid, Obj.desg)

        if isinstance(obj, str):
            query = query.filter(Obj.desg == obj)
        else:
            query = query.filter(Obj.objid == obj)

        try:
            objid, desg = query.one()
        except NoResultFound:
            if isinstance(obj, str):
                objid = None
                desg = obj
            else:
                raise BadObjectID('{} not found in database'.format(obj))

        return objid, desg
