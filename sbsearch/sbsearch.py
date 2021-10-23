# Licensed with the 3-clause BSD license.  See LICENSE for details.

__all__ = ['SBSearch']

from typing import Any, Dict, List, Optional, Tuple, TypeVar, Union
from logging import Logger

import numpy as np
import sqlalchemy as sa
from sqlalchemy.orm import Session, Query
from sqlalchemy.engine import reflection
from astropy.time import Time

from .ephemeris import get_ephemeris_generator, EphemerisGenerator
from .sbsdb import SBSDatabase
from .model import (Base, Ephemeris, Observation, Found)
from .spatial import (  # pylint: disable=E0611
    SpatialIndexer, polygon_string_intersects_line,
    polygon_string_intersects_about_line)
from .target import MovingTarget
from .exceptions import DesignationError, UnknownSource
from .config import Config
from .logging import setup_logger, ProgressBar


SBSearchObject = TypeVar('SBSearchObject', bound='SBSearch')


class SBSearch:
    """Small-body search tool.


    Parameters
    ----------
    database : string or sqlalchemy Session
        The sqlalchemy-formatted database URL or a sqlalchemy session
        to use.

    min_edge_length : float, optional
        Minimum edge length to index, radians.  See
        http://s2geometry.io/resources/s2cell_statistics for cell sizes.

    uncertainty_ellipse : bool, optional
        Search considering the uncertainty ellipse.

    padding : float, optional
        Additional padding to the search area, arcmin.

    *args
        Optional `SBSDatabase` arguments.

    """

    def __init__(self, database: Union[str, Session], *args,
                 min_edge_length: float = 3e-4, padding: float = 0,
                 uncertainty_ellipse: bool = False,
                 log: str = '/dev/null', logger_name: str = 'SBSearch'
                 ) -> None:
        self.db = SBSDatabase(database, *args)
        self.db.verify()
        self.indexer: SpatialIndexer = SpatialIndexer(min_edge_length)
        self._source: Union[Observation, None] = None
        self.uncertainty_ellipse: bool = uncertainty_ellipse
        self.padding: float = padding
        self.logger: Logger = setup_logger(filename=log, name=logger_name)

    def __enter__(self) -> SBSearchObject:
        return self

    def __exit__(self, *args):
        self.db.session.commit()
        self.db.session.close()
        self.logger.info('Terminated')

    @classmethod
    def with_config(cls, config: Config) -> SBSearchObject:
        """Instantiate with these configuration options."""
        return cls(**config.config)

    @property
    def source(self) -> Observation:
        """Observation data source for searches.

        Only this survey data source will be searched.

        Raises
        ------
        ValueError if the source has not been set.

        """
        if self._source is None:
            raise ValueError('sbsearch data source not set.')

        return self._source

    @source.setter
    def source(self, source: Union[str, Observation]) -> None:
        """Set the observation data source for searches.

        May be set to a table name, or a data model object, derived from
        ``model.Observation``, e.g.,

        >> > from sbsearch.model import ExampleSurvey
        >> > sbs.source = ExampleSurvey
        >> > print(ExampleSurvey)
        'example_survey'
        >> > sbs.source = 'example_survey'

        Or, to search all data, regardless of source:
        >> > from sbsearch.model import Observation
        >> > sbs.source = Observation

        But note that in this case moving target ephemerides will not change
        with observatory location, using ``Observation.__obscode__`` for all
        data sources.

        """

        if isinstance(source, str):
            e: Exception
            try:
                self._source = self.sources[source]
            except KeyError as e:
                raise UnknownSource(source) from e
        else:
            if source not in self.sources.values():
                raise UnknownSource(source)
            self._source = source

    @property
    def sources(self) -> Dict[str, Observation]:
        """Dictionary of observation data sources in the information model.

        The dictionary is keyed by database table name.

        """
        return {source.__tablename__: source
                for source in Observation.__subclasses__()}

    @property
    def uncertainty_ellipse(self) -> bool:
        """Set to search the uncertainty ellipse.

        This is in addition to ``padding``.

        """
        return self._uncertainty_ellipse

    @uncertainty_ellipse.setter
    def uncertainty_ellipse(self, flag: bool):
        self._uncertainty_ellipse: bool = flag

    @property
    def padding(self) -> float:
        """Set to pad ephemeris regions by this amount in arcmin.

        This is in addition to ``uncertainty_ellipse``.

        """
        return self._padding

    @padding.setter
    def padding(self, amount: float):
        self._padding: float = max(amount, 0)

    def _get_spatial_term_class(self, obs: Union[Observation, Base, None] = None) -> Base:
        if obs is None:
            return self.source.terms.property.mapper.class_
        elif isinstance(obs, type):
            return obs.terms.property.mapper.class_
        else:
            return type(obs).terms.property.mapper.class_

    def _get_spatial_term_index_name(self, SpatialTerm: Base) -> str:
        inspector: reflection.Inspector = reflection.Inspector.from_engine(
            self.db.engine)
        SpatialTermIndex: sa.Index
        for index in inspector.get_indexes(SpatialTerm.__tablename__):
            if index['column_names'] == ['term']:
                return index['name']

    def add_designation(self, designation: str) -> MovingTarget:
        """Add designation to database and return moving target.


        Parameters
        ----------
        designation : str
            The target's designation.


        Returns
        -------
        target: MovingTarget
            The newly created target.

        """

        target: MovingTarget = MovingTarget(designation, db=self.db)
        target.add()
        return target

    def get_designation(self, designation: str, add: bool = False) -> MovingTarget:
        """Get target named ``designation`` from database.


        Parameters
        ----------
        designation : str
            The target designation.

        add : bool, optional
            If the target does not exist, add it.

        """
        try:
            return MovingTarget.from_designation(designation, db=self.db)
        except DesignationError:
            if add:
                return self.add_designation(designation)
            raise

    def get_object_id(self, object_id: int) -> MovingTarget:
        """Get target by database object ID."""
        return MovingTarget.from_id(object_id, db=self.db)

    def add_ephemeris(self, observer: str, target: MovingTarget,
                      start_date: str, stop_date: str, cache: bool = False
                      ) -> None:
        """Add ephemeris to database.

        Parameters
        ----------
        observer: str
            Observatory code or location.  Must be resolvable by the
            ephemeris generator.

        target: MovingTarget
            Must have an object ID in the database and be resolvable by the
            ephemeris generator.

        start_date, stop_date: str
            Start and stop dates, UTC, in a format parseable by astropy `Time`.

        cache: bool, optional
            Use cached results, if possible, otherwise cache results.  For
            ephemerides generated via astroquery.

        """

        start = Time(start_date)
        stop = Time(stop_date)
        g: EphemerisGenerator = get_ephemeris_generator()
        ephemerides: List[Ephemeris] = g.target_over_date_range(
            observer, target, start, stop, cache=cache)
        eph: Ephemeris
        for eph in ephemerides:
            self.db.session.add(eph)
        self.db.session.commit()
        self.logger.info('Added %d ephemeris point%s for %s at %s.',
                         len(ephemerides),
                         '' if len(ephemerides) == 1 else 's',
                         target.primary_designation, observer)

    def get_ephemeris(self, target: MovingTarget, start_date: str,
                      stop_date: str) -> List[Ephemeris]:
        """Get ephemeris from database.

        Parameters
        ----------
        target: MovingTarget
            Must have an object ID in the database and be resolvable by the
            ephemeris generator.

        start_date, stop_date: str
            Start and stop dates, UTC, in a format parseable by astropy `Time`.

        Returns
        -------
        eph: list of Ephemeris objects

        """

        start = Time(start_date)
        stop = Time(stop_date)

        # expand date search by small amount to facilitate floating point
        # comparisons
        return (
            self.db.session.query(Ephemeris)
            .filter(Ephemeris.object_id == target.object_id)
            .filter(Ephemeris.mjd >= (start.mjd - 0.000001))
            .filter(Ephemeris.mjd <= (stop.mjd + 0.000001))
            .all()
        )

    def add_observations(self, observations: List[Observation]):
        """Add observations to the database.

        Observations must be added with this method, or else the spatial index
        will not include them.

        """

        for obs in observations:
            self.db.session.add(obs)
            # get spatial term class specific to this observation class
            SpatialTerm: Any = self._get_spatial_term_class(obs)
            for term in self.indexer.index_polygon_string(obs.fov):
                obs.terms.append(
                    SpatialTerm(
                        source_id=obs.id,
                        term=term
                    ))

        self.db.session.commit()
        self.logger.info('Added %d observation%s.', len(observations),
                         '' if len(observations) == 1 else 's')

    def get_observations(self, source: Optional[str] = None,
                         mjd: Optional[List[float]] = None
                         ) -> List[Observation]:
        """Get observations from database.

        Parameters
        ----------
        source: str, optional
            Get observations from this source.

        mjd: list of float, optional
            Get observations between these modified Julian dates.

        """

        q: Query = self.db.session.query(self.source)
        if source is not None:
            q = q.filter(Observation.source == source)
        if mjd is not None:
            q = q.filter(Observation.mjd_start <= max(mjd)).filter(
                Observation.mjd_stop >= min(mjd))

        return q.all()

    def re_index(self, drop_index: bool = True):
        """Delete and recreate the spatial index for the current source.

        To change the minimum edge length of a database, first initialize
        SBSearch with the new value, then call this method for all
        observations.


        Parameters
        ----------
        drop_index : bool, optional
            Speed up inserts by dropping the database index, and rebuilding
            it after adding the new terms.

        """

        SpatialTerm: Base = self._get_spatial_term_class(self.source)
        spatialtermindex_name: sa.Index = self._get_spatial_term_index_name(
            SpatialTerm)

        n_terms: int = self.db.session.query(SpatialTerm).delete(
            synchronize_session=False
        )
        self.db.session.commit()
        self.logger.info('Deleted %d spatial term%s.', n_terms,
                         '' if n_terms == 1 else 's')

        if drop_index:
            # drop the database index to speed up adding many terms
            self.db.session.execute(f'DROP INDEX {spatialtermindex_name}')

        n_terms = 0
        n_obs: int = self.db.session.query(self.source).count()
        with ProgressBar(n_obs, self.logger, scale='log') as bar:
            n_obs = 0
            while True:
                observations: Query = (
                    self.db.session.query(self.source)
                    .offset(n_obs)
                    .limit(1000)
                    .all()
                )
                if len(observations) == 0:
                    break

                for obs in observations:
                    n_obs += 1
                    bar.update()
                    for term in self.indexer.index_polygon_string(obs.fov):
                        n_terms += 1
                        self.db.session.add(
                            SpatialTerm(
                                source_id=obs.id,
                                term=term
                            ))

                # check-in to avoid soaking up too much memory
                self.db.session.commit()

        self.db.session.commit()

        if drop_index:
            self.logger.info('Rebuilding database index.')
            self.db.session.execute(
                f'CREATE INDEX {spatialtermindex_name} '
                f'ON {spatialtermindex_name[3:]} (term)'
            )

        self.logger.info('Re-indexed %d observation%s with %d spatial term%s.',
                         n_obs, '' if n_obs == 1 else 's',
                         n_terms, '' if n_terms == 1 else 's')

    def add_found(self, target: MovingTarget, observations: List[Observation],
                  cache: bool = True) -> None:
        """Add observations of a target to the found database.


        Parameters
        ----------
        target : MovingTarget
            The found target.  Must already exist in the database.

        observations : list of Observation
            The observations in which the target is found.  Must all be
            from the same source.

        cache: bool, optional
            Use cached results, if possible, otherwise cache results.  For
            ephemerides generated via astroquery.


        Returns
        -------
        found : list of Found
            The inserted items.

        """

        # get moving target metadata
        g: EphemerisGenerator = get_ephemeris_generator()

        # verify that all sources are the same
        if len(set([obs.source for obs in observations])) > 1:
            raise ValueError('all observations must be from the same source')

        found: List[Found] = []
        observer = observations[0].__obscode__
        dates: Time = Time([(obs.mjd_start + obs.mjd_stop) / 2
                            for obs in observations], format='mjd')
        ephemerides: List[Ephemeris] = g.target_at_dates(
            observer, target, dates, cache=cache)
        for eph, obs in zip(ephemerides, observations):
            f: Found = Found(
                object_id=target.object_id,
                observation_id=obs.observation_id,
            )
            for k in ['mjd', 'rh', 'delta', 'phase', 'drh', 'true_anomaly',
                      'ra', 'dec', 'dra', 'ddec', 'unc_a', 'unc_b',
                      'unc_theta', 'elong', 'sangle', 'vangle', 'vmag',
                      'retrieved']:
                setattr(f, k, getattr(eph, k))

            self.db.session.add(f)
            found.append(f)

        self.db.session.commit()
        self.logger.info('Added %d observation%s of %s.', len(found),
                         '' if len(found) == 1 else 's',
                         target.primary_designation)
        return found

    def get_found(self, target: Optional[MovingTarget] = None,
                  mjd: Optional[List[float]] = None
                  ) -> List[Any]:
        """Get found objects from database.


        Parameters
        ----------
        target : MovingTarget, optional
            Get found observations of this target.

        mjd : list of float
            Get found observations between these modified Julian dates.

        join : list of table objects
            Join rows with these table objects, e.g., 'Obj', 'Observation'.

        """

        q: Query = self.db.session.query(Found)

        if target is not None:
            q = q.filter(Found.object_id == target.object_id)

        if mjd is not None:
            q = (q.filter(Found.mjd >= min(mjd))
                 .filter(Found.mjd <= max(mjd)))

        return q.all()

    def find_observations_intersecting_polygon(
        self, ra: np.ndarray, dec: np.ndarray
    ) -> List[Observation]:
        """Find observations intersecting given polygon.

        Parameters
        ----------
        vertices: array-like
            Polygon vertices in units of radians.  It is assumed that the
            area is < 2 pi steradians.

        Returns
        -------
        observations: list of Observation

        """

        terms: List[str] = self.indexer.query_polygon(
            np.array(ra, float), np.array(dec, float))
        SpatialTerm: Base = self._get_spatial_term_class()
        obs: List[Observation] = (
            self.db.session.query(self.source)
            .join(SpatialTerm)
            .filter(SpatialTerm.term.in_(terms))
            .all()
        )
        return obs

    def find_observations_intersecting_line(
        self, ra: np.ndarray, dec: np.ndarray,
        a: Optional[np.ndarray] = None, b: Optional[np.ndarray] = None,
        approximate: bool = False
    ) -> List[Observation]:
        """Find observations intersecting given line.

        Parameters
        ----------
        ra, dec: array-like
            Line vertices in units of radians.

        a, b: array-like, optional
            Extend the search area about the line by these angular distances,
            radians.  ``a`` is parallel to the line segment, ``b`` is
            perpendicular.  Only the first and last ``a`` are considered.

        approximate: bool, optional
            Do not check potential matches in detail.

        Returns
        -------
        observations: list of Observation

        """

        # normalize inputs for use with spatial submodule
        _ra = np.array(ra, float)
        _dec = np.array(dec, float)
        _a: Optional[np.ndarray] = None
        _b: Optional[np.ndarray] = None
        query_about: bool = False
        if a is not None and b is not None:
            _a: np.ndarray = np.array(a, float)
            _b: np.ndarray = np.array(b, float)
            if len(_b) != len(_a):
                raise ValueError('Size of a and be must match.')
            query_about = True

        obs: List[Observation] = []
        # search for the full line at once
        terms: List[str]
        if query_about:
            terms = self.indexer.query_about_line(
                _ra, _dec, _a, _b)[0]
        else:
            terms = self.indexer.query_line(_ra, _dec)

        SpatialTerm: Base = self._get_spatial_term_class()
        _obs: List[Observation] = (
            self.db.session.query(self.source)
            .join(SpatialTerm)
            .filter(SpatialTerm.term.in_(terms))
            .all()
        )

        if not approximate:
            # test each observation for intersection with the observation FOV
            for o in _obs:
                if query_about:
                    intersects = polygon_string_intersects_about_line(
                        o.fov, _ra, _dec, _a, _b)
                else:
                    intersects = polygon_string_intersects_line(
                        o.fov, _ra, _dec)
                if intersects:
                    obs.append(o)
        else:
            obs = _obs

        return obs

    @ staticmethod
    def _test_line_intersection_with_observations_at_time(
        obs: List[Observation], ra: np.ndarray, dec: np.ndarray,
        mjd: np.ndarray, a: Optional[np.ndarray] = None,
        b: Optional[np.ndarray] = None
    ) -> List[Observation]:
        """Test each observation for intersection with the observation FOV."""
        query_about: bool = a is not None
        _obs: List[Observation] = []
        for o in obs:
            dt = mjd[1] - mjd[0]
            line_start = (o.mjd_start - mjd[0]) / dt
            line_stop = (o.mjd_stop - mjd[0]) / dt
            if query_about:
                intersects = polygon_string_intersects_about_line(
                    o.fov, ra, dec, a, b,
                    line_start=line_start, line_stop=line_stop
                )
            else:
                intersects = polygon_string_intersects_line(
                    o.fov, ra, dec,
                    line_start=line_start, line_stop=line_stop)
            if intersects:
                _obs.append(o)

        return _obs

    def find_observations_intersecting_line_at_time(
        self, ra: np.ndarray, dec: np.ndarray, mjd: np.ndarray,
        a: Optional[np.ndarray] = None, b: Optional[np.ndarray] = None,
        approximate: bool = False
    ) -> List[Observation]:
        """Find observations intersecting given line at given times.

        Each segment is separately queried.

        Parameters
        ----------
        ra, dec: array-like
            Line vertices in units of radians.

        mjd: array-like
            Search by time for each point, UTC.

        a, b: array-like, optional
            Extend the search area about the line by these angular distances,
            radians.  ``a`` is parallel to the line segment, ``b`` is
            perpendicular.

        approximate: bool, optional
            Do not check potential matches in detail.

        Returns
        -------
        observations: list of Observation

        """

        # normalize inputs for use with spatial submodule
        _ra = np.array(ra, float)
        _dec = np.array(dec, float)

        N: int = len(_ra)
        if len(_dec) != N:
            raise ValueError('ra and dec must have same length')

        _a: Optional[np.ndarray] = None
        _b: Optional[np.ndarray] = None
        if a is not None or b is not None:
            _a: np.ndarray = np.array(a, float)
            _b: np.ndarray = np.array(b, float)
            if len(_a) != N or len(_b) != N:
                raise ValueError('ra, dec, a, and b must have same length')

        SpatialTerm: Base = self._get_spatial_term_class()
        obs: List[Observation] = []
        i: int
        for i in range(N - 1):
            query_about: bool = a is not None

            terms: List[str]
            if query_about:
                terms = self.indexer.query_about_line(
                    _ra[i:i + 2], _dec[i:i + 2], _a[i:i + 2], _b[i:i + 2])[0]
            else:
                terms = self.indexer.query_line(_ra[i:i + 2], _dec[i:i + 2])

            nearby_obs: List[Observation] = (
                self.db.session.query(self.source)
                .join(SpatialTerm)
                .filter(SpatialTerm.term.in_(terms))
                .filter(Observation.mjd_start <= mjd[i + 1])
                .filter(Observation.mjd_stop >= mjd[i])
                .all()
            )

            if len(nearby_obs) == 0:
                continue
            elif approximate:
                obs.extend(nearby_obs)
            else:
                # check for detailed intersection
                obs.extend(
                    self._test_line_intersection_with_observations_at_time(
                        nearby_obs, _ra[i:i + 2], _dec[i:i + 2], mjd[i:i + 2],
                        None if a is None else _a[i:i + 2],
                        None if b is None else _b[i:i + 2]
                    )
                )

        # duplicates can accumulate because each segment is searched
        # individually
        return list(set(obs))

    def find_observations_by_ephemeris(self, eph: List[Ephemeris]
                                       ) -> List[Observation]:
        """Find observations covering given ephemeris.


        Parameters
        ----------
        eph: list of Ephemeris
            The ephemeris points to check.  Assumed to be continuous.


        Returns
        -------
        obs: list of Observation


        Notes
        ------
        See attributes ``uncertainty_ellipse`` and ``padding`` for search
        options.

        If searching over the uncertainty ellipse, the uncertainty area is
        circumscribed with a quadrilateral.

        """

        N: int = len(eph)

        # unpack ephemerides into arrays
        ra: np.ndarray
        dec: np.ndarray
        mjd: np.ndarray
        unc_a: np.ndarray
        unc_b: np.ndarray
        unc_theta: np.ndarray
        ra, dec, mjd, unc_a, unc_b, unc_theta = np.empty((6, N))
        for i in range(N):
            ra[i] = eph[i].ra  # deg
            dec[i] = eph[i].dec  # deg
            mjd[i] = eph[i].mjd  # UTC
            # enforce a minimum uncertainty to avoid overlapping vertices
            unc_a[i] = 1 if eph[i].unc_a is None else eph[i].unc_a  # arcsec
            unc_b[i] = 1 if eph[i].unc_b is None else eph[i].unc_b  # arcsec
            # deg:
            unc_theta[i] = 0 if eph[i].unc_theta is None else eph[i].unc_theta

        # convert to radians
        ra = np.radians(ra)
        dec = np.radians(dec)
        unc_a = np.radians(unc_a / 3600)
        unc_b = np.radians(unc_b / 3600)
        padding: np.ndarray = np.ones_like(ra) * np.radians(self.padding / 60)

        obs: List[Observation]
        if self.uncertainty_ellipse:
            a: np.ndarray
            b: np.ndarray
            a, b = self._ephemeris_uncertainty_offsets(eph)
            obs = self.find_observations_intersecting_line_at_time(
                ra, dec, mjd, a=a + padding, b=b + padding
            )
        elif self.padding > 0:
            obs = self.find_observations_intersecting_line_at_time(
                ra, dec, mjd, a=padding, b=padding
            )
        else:
            obs = self.find_observations_intersecting_line_at_time(
                ra, dec, mjd
            )

        self.logger.info('%d observation%s found.', len(obs),
                         '' if len(obs) == 1 else 's')

        return obs

    @ staticmethod
    def _ephemeris_uncertainty_offsets(eph: List[Ephemeris]
                                       ) -> Tuple[np.ndarray]:
        """Generate ephemeris offsets that cover the uncertainty area.

        Requires the following definitions in the Ephemeris object:
            dra, ddec, unc_a, unc_b, unc_theta


        Parameters
        ----------
        eph: list of Ephemeris
            Must be at least 2 points.


        Returns
        -------
        a, b: np.ndarray
            The offsets, suitable for ``find_observations_intersecting_line``.

        """

        if len(eph) < 2:
            raise ValueError('Must have at least 2 ephemeris points.')

        dra: np.ndarray = np.radians([e.dra for e in eph])
        ddec: np.ndarray = np.radians([e.ddec for e in eph])
        unc_a: np.ndarray = np.radians([e.unc_a / 3600 for e in eph])
        unc_b: np.ndarray = np.radians([e.unc_b / 3600 for e in eph])
        unc_theta: np.ndarray = np.radians([e.unc_theta for e in eph])

        # target motion as position angle, radians
        pa: np.ndarray = np.arctan2(dra, ddec)

        # set up proper motion unit vectors: R in the direction of motion,
        # P perpendicular to it
        R: np.ndarray = np.array((np.cos(pa), np.sin(pa)))
        P: np.ndarray = np.array((np.cos(pa + np.pi / 2),
                                  np.sin(pa + np.pi / 2)))

        # setup uncertainty vectors
        A: np.ndarray = unc_a * np.array((np.cos(unc_theta),
                                          np.sin(unc_theta)))
        B: np.ndarray = unc_b * np.array((np.cos(unc_theta + np.pi / 2),
                                          np.sin(unc_theta + np.pi / 2)))

        # compute offsets a, b
        a = np.max(np.abs((np.sum(A * R, 0), np.sum(B * R, 0))), 0)
        b = np.max(np.abs((np.sum(A * P, 0), np.sum(B * P, 0))), 0)

        return a, b
