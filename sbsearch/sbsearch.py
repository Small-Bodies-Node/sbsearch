# Licensed with the 3-clause BSD license.  See LICENSE for details.

__all__ = ['SBSearch']

from typing import Any, Dict, List, Optional, Tuple, TypeVar, Union
from logging import Logger, getLogger

import numpy as np
from sqlalchemy.orm import Session
from astropy.time import Time

from .ephemeris import get_ephemeris_generator, EphemerisGenerator
from .sbsdb import SBSDatabase
from .model import Ephemeris, Observation, ObservationSpatialTerm
from .spatial import (SpatialIndexer, polygon_string_intersects_line,
                      polygon_string_intersects_about_line)
from .target import MovingTarget
from .exceptions import UnknownSource
from .config import Config
from .logging import setup_logger


SBSearchObject = TypeVar('SBSearchObject', bound='SBSearch')


class SBSearch:
    """Small-body search tool.


    Parameters
    ----------
    min_edge_length : float
        Minimum edge length to index, radians.  See
        http://s2geometry.io/resources/s2cell_statistics for cell sizes.


    database : string or sqlalchemy Session
        The sqlalchemy-formatted database URL or a sqlalchemy session
        to use.

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
        self.source: Observation = Observation
        self.uncertainty_ellipse: bool = uncertainty_ellipse
        self.padding: float = padding
        self.logger: Logger = setup_logger(filename=log, name=logger_name)

    def __enter__(self) -> SBSearchObject:
        return self

    def __exit__(self, *args):
        self.db.session.commit()
        self.db.session.close()
        self.logger.info('Terminated at %sZ', Time.now().iso)

    @ classmethod
    def with_config(cls, config: Config) -> SBSearchObject:
        """Instantiate with these configuration options."""
        return cls(**config.config)

    @ property
    def source(self) -> Observation:
        """Observation data source for searches.

        Only this survey data source will be searched.

        """
        return self._source

    @ source.setter
    def source(self, source: Union[str, Observation]) -> None:
        """Set the observation data source for searches.

        May be set to a table name, or a data model object, derived from
        ``model.Observation``, e.g.,

        >> > from sbsearch.model import UnspecifiedSurvey
        >> > sbs.source = UnspecifiedSurvey
        >> > print(UnspecifiedSurvey)
        'unspecified_survey'
        >> > sbs.source = 'unspecified_survey'

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

    @ property
    def sources(self) -> Dict[str, Observation]:
        """Dictionary of observation data sources in the information model.

        The dictionary is keyed by database table name.

        """
        return {source.__tablename__: source
                for source in [Observation] + Observation.__subclasses__()}

    @ property
    def uncertainty_ellipse(self) -> bool:
        """Set to search the uncertainty ellipse.

        This is in addition to ``padding``.

        """
        return self._uncertainty_ellipse

    @ uncertainty_ellipse.setter
    def uncertainty_ellipse(self, flag: bool):
        self._uncertainty_ellipse: bool = flag

    @ property
    def padding(self) -> float:
        """Set to pad ephemeris regions by this amount in arcmin.

        This is in addition to ``uncertainty_ellipse``.

        """
        return self._padding

    @ padding.setter
    def padding(self, amount: float):
        self._padding: float = amount

    def add_designation(self, designation: str) -> MovingTarget:
        """Add designation to database and return moving target."""
        target: MovingTarget = MovingTarget(designation, db=self.db)
        target.add()
        return target

    def get_designation(self, designation: str) -> MovingTarget:
        """Get target named ``designation`` from database."""
        return MovingTarget.from_designation(designation, db=self.db)

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
        eph: Ephemeris
        for eph in g.target_over_date_range(observer, target, start, stop,
                                            cache=cache):
            self.db.session.add(eph)
        self.db.session.commit()

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
            for term in self.indexer.index_polygon_string(obs.fov):
                obs.terms.append(
                    ObservationSpatialTerm(
                        observation_id=obs.observation_id,
                        term=term
                    ))

        self.db.session.commit()

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

        q: Any = self.db.session.query(Observation)
        if source is not None:
            q = q.filter(Observation.source == source)
        if mjd is not None:
            q = q.filter(Observation.mjd_start <= max(mjd)).filter(
                Observation.mjd_stop >= min(mjd))

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
        obs: List[Observation] = (
            self.db.session.query(Observation)
            .join(ObservationSpatialTerm)
            .filter(ObservationSpatialTerm.term.in_(terms))
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
        if a is not None or b is not None:
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

        _obs: List[Observation] = (
            self.db.session.query(Observation)
            .join(ObservationSpatialTerm)
            .filter(ObservationSpatialTerm.term.in_(terms))
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
        approximate: bool = False,
        source: Optional[Union[str, Observation]] = None
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
                self.db.session.query(Observation)
                .join(ObservationSpatialTerm)
                .filter(ObservationSpatialTerm.term.in_(terms))
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
        See ``uncertainty_ellipse`` and ``padding`` for search options.

        If searching over the uncertainty ellipse, the area is approximated by
        a polygonal region.

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
