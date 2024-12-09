# Licensed with the 3-clause BSD license.  See LICENSE for details.

__all__ = ["SBSearch", "IntersectionType"]

import enum
from typing import Any, Dict, List, Optional, Tuple, TypeVar, Union
import logging
from logging import Logger

import numpy as np
from sqlalchemy.orm import Session, Query
from astropy.time import Time
from astropy.coordinates import SkyCoord

from . import core
from .ephemeris import get_ephemeris_generator, EphemerisGenerator
from .sbsdb import SBSDatabase
from .model import Base, Ephemeris, Observation, Found  # noqa F203
from .spatial import (  # pylint: disable=E0611
    SpatialIndexer,
    polygon_intersects_line,
    polygon_intersects_about_line,
    polygon_intersects_polygon,
    polygon_contains_point,
    polygon_intersects_cap,
)
from .target import MovingTarget, FixedTarget
from .exceptions import DesignationError, UnknownSource
from .config import Config
from .logging import ProgressTriangle, setup_logger, setup_search_logger


# keep synced with libspatial.cpp, test_spatial.py
class IntersectionType(enum.Enum):
    ImageContainsCenter = 0
    ImageContainsArea = 1
    ImageIntersectsArea = 2
    AreaContainsImage = 3


SBSearchObject = TypeVar("SBSearchObject", bound="SBSearch")


class SBSearch:
    """Small-body search tool.

    Two loggers are used.  The first, under ``logger_name``, publishes general
    messages to the console and to a log file (``log``).  The second,
    ``logger_name (search)`` publishes moving target search progress to the
    console.


    Parameters
    ----------
    database : string or sqlalchemy Session
        The sqlalchemy-formatted database URL or a sqlalchemy session to use.

    min_edge_length : float, optional max_edge_length : float, optional
        Minimum and maximum edge length to index, radians.  See
        http://s2geometry.io/resources/s2cell_statistics for cell sizes (1
        radian = 6380 km on the Earth).

    uncertainty_ellipse : bool, optional
        Search considering the uncertainty ellipse.

    padding : float, optional
        Additional padding to the search area, arcmin.

    log : str, optional
        Log file name.

    logger_name : str, optional
        Use this logger instance name.

    debug : bool, optional
        Enable debugging messages.

    *args
        Optional `SBSDatabase` arguments.

    """

    def __init__(
        self,
        database: Union[str, Session],
        *args,
        min_edge_length: float = 3e-4,
        max_edge_length: float = 0.017,
        uncertainty_ellipse: bool = False,
        padding: float = 0,
        start_date: Optional[Time] = None,
        stop_date: Optional[Time] = None,
        intersection_type: IntersectionType = IntersectionType.ImageIntersectsArea,
        log: str = "/dev/null",
        logger_name: str = "SBSearch",
        arc_limit: float = 0.17,
        time_limit: float = 365,
        debug: bool = False,
    ) -> None:
        self.db = SBSDatabase(database, *args)
        self.db.verify()
        self.indexer: SpatialIndexer = SpatialIndexer(min_edge_length, max_edge_length)
        self._source: Union[Observation, None] = None
        self.uncertainty_ellipse: bool = uncertainty_ellipse
        self.padding: float = padding
        self.start_date: Union[Time, None] = start_date
        self.stop_date: Union[Time, None] = stop_date
        self.intersection_type: IntersectionType = intersection_type
        self.arc_limit = arc_limit
        self.time_limit = time_limit
        self.debug = debug
        log_level: int = logging.DEBUG if debug else None
        self.logger: Logger = setup_logger(
            filename=log, name=logger_name, level=log_level
        )
        self.search_logger: Logger = setup_search_logger(logger_name)

    def __enter__(self) -> SBSearchObject:
        return self

    def __exit__(self, *args):
        self.db.session.commit()
        self.db.session.close()
        self.logger.info("Terminated")

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
            raise ValueError("sbsearch data source not set.")

        return self._source

    @source.setter
    def source(self, source: Union[str, Observation]) -> None:
        """Set the observation data source for searches.

        May be set to a table name, or a data model object, derived from
        ``model.Observation``, e.g.,

            from sbsearch.model import ExampleSurvey
            sbs.source = ExampleSurvey
            print(ExampleSurvey)
            -> 'example_survey'
            sbs.source = 'example_survey'

        Or, to search all data, regardless of source:

            from sbsearch.model import Observation
            sbs.source = Observation

        But note that in this case moving target ephemerides will not change
        with observatory location, using ``Observation.__obscode__`` for all
        data sources.

        """

        if isinstance(source, str):
            if source == "observation":
                self._source = Observation
            else:
                e: Exception
                try:
                    self._source = self.sources[source]
                except KeyError as e:  # noqa F841
                    raise UnknownSource(source) from e
        else:
            if source == Observation:
                pass
            elif source not in self.sources.values():
                raise UnknownSource(source)
            self._source = source

    @property
    def sources(self) -> Dict[str, Observation]:
        """Dictionary of observation data sources in the information model.

        The dictionary is keyed by database table name.

        """
        return {source.__tablename__: source for source in Observation.__subclasses__()}

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

    def add_ephemeris(
        self,
        observer: str,
        target: MovingTarget,
        start_date: str,
        stop_date: str,
        cache: bool = False,
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
            observer, target, start, stop, cache=cache
        )
        eph: Ephemeris
        for eph in ephemerides:
            self.db.session.add(eph)
        self.db.session.commit()
        self.logger.info(
            "Added %d ephemeris point%s for %s at %s.",
            len(ephemerides),
            "" if len(ephemerides) == 1 else "s",
            target.primary_designation,
            observer,
        )

    def get_ephemeris(
        self, target: MovingTarget, start_date: str, stop_date: str
    ) -> List[Ephemeris]:
        """Get ephemeris from database.

        Parameters
        ----------
        target: MovingTarget
            Must have an object ID in the database and be resolvable by the
            ephemeris generator.

        start_date, stop_date: str
            Start and stop dates, UTC, in a format parsable by astropy `Time`.

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

    def add_observations(self, observations: List[Observation]) -> None:
        """Add observations to the database.

        If ``spatial_terms`` is not set, then new terms are generated.


        Parameters
        ----------
        observations: list of Observations

        """

        for obs in observations:
            if obs.spatial_terms is None:
                obs.spatial_terms = self.indexer.index_polygon(
                    *core.polygon_string_to_arrays(obs.fov)
                )

        self.db.session.add_all(observations)
        self.db.session.commit()

        self.logger.debug(
            "Added %d observation%s.",
            len(observations),
            "" if len(observations) == 1 else "s",
        )

    def update_observations(
        self, observations: List[Observation], reindex: bool | None = True
    ) -> None:
        """Update observations in the database.

        Parameters
        ----------
        observations: list of Observations
            The observations to update.  They must have observation_id defined.

        reindex: bool, optional
            If ``True`` then the spatial index terms will be regenerated.

        """

        if reindex:
            for obs in observations:
                obs.spatial_terms = self.indexer.index_polygon(
                    *core.polygon_string_to_arrays(obs.fov)
                )

        for obs in observations:
            self.db.session.merge(obs)
        self.db.session.commit()

        self.logger.debug(
            "Updated %d observation%s.",
            len(observations),
            "" if len(observations) == 1 else "s",
        )

    def get_observations(self) -> List[Observation]:
        """Get observations from database."""
        q: Query = self.db.session.query(self.source)
        q = self._filter_by_source(q)
        q = self._filter_by_date(q)
        return q.all()

    def re_index(self, terms=True) -> None:
        """Delete and recreate the spatial index for the current source.

        To change the minimum edge length of a database, first initialize
        SBSearch with the new value, set the data source to `Observation`,
        then call this method with ``terms=True``.


        Parameters
        ----------
        terms : bool, optional
            If ``False``, do not regenerate the spatial index terms, only
            regenerate the database index of the terms.

        """

        n_obs: int = self.db.session.query(self.source).count()
        if terms:
            self.logger.info(
                "Generating spatial index terms for %d %s rows.",
                n_obs,
                self.source.__tablename__,
            )
        else:
            self.logger.info("Recreating spatial term index.")

        self.db.drop_spatial_index()

        if terms:
            n_terms: int = 0
            with ProgressTriangle(1, self.logger, base=2) as tri:
                n_obs = 0
                while True:
                    obs: Observation
                    count: int = 0
                    for obs in (
                        self.db.session.query(self.source)
                        .order_by(self.source.observation_id)
                        .offset(n_obs)
                        .limit(10000)
                    ):
                        count += 1
                        n_obs += 1
                        tri.update()
                        terms = self.indexer.index_polygon(
                            *core.polygon_string_to_arrays(obs.fov)
                        )
                        n_terms += len(terms)
                        obs.spatial_terms = terms

                    if count == 0:
                        break

                    # check-in to avoid soaking up too much memory
                    self.db.session.commit()

        self.db.create_spatial_index()
        self.db.session.commit()

        if terms:
            self.logger.info(
                "Re-indexed %d observation%s with %d spatial term%s.",
                n_obs,
                "" if n_obs == 1 else "s",
                n_terms,
                "" if n_terms == 1 else "s",
            )
        else:
            self.logger.info("Indexing complete.")

    def add_found(
        self,
        target: MovingTarget,
        observations: List[Observation],
        cache: bool = True,
    ) -> None:
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
            raise ValueError("all observations must be from the same source")

        found: List[Found] = []
        observer = observations[0].__obscode__
        dates: Time = Time(
            [(obs.mjd_start + obs.mjd_stop) / 2 for obs in observations],
            format="mjd",
        )
        ephemerides: List[Ephemeris] = g.target_at_dates(
            observer, target, dates, cache=cache
        )
        for eph, obs in zip(ephemerides, observations):
            f: Found = Found(
                object_id=target.object_id,
                observation_id=obs.observation_id,
            )
            for k in [
                "mjd",
                "rh",
                "delta",
                "phase",
                "drh",
                "true_anomaly",
                "ra",
                "dec",
                "dra",
                "ddec",
                "unc_a",
                "unc_b",
                "unc_theta",
                "elong",
                "sangle",
                "vangle",
                "vmag",
                "retrieved",
            ]:
                setattr(f, k, getattr(eph, k))

            self.db.session.add(f)
            found.append(f)

        self.db.session.commit()
        self.logger.info(
            "Added %d observation%s of %s.",
            len(found),
            "" if len(found) == 1 else "s",
            target.primary_designation,
        )
        return found

    def get_found(
        self,
        target: Optional[MovingTarget] = None,
        mjd: Optional[List[float]] = None,
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
            q = q.filter(Found.mjd >= min(mjd)).filter(Found.mjd <= max(mjd))

        return q.all()

    def _filter_by_source(self, query):
        """Filter query by observation source."""

        if self.source != Observation:
            query = query.filter(Observation.source == self.source.__tablename__)
        else:
            # If the search is for any observation, we need to limit the results
            # to sources that are known to this version of sbsearch as the
            # database, especially the dev database, may have surveys unknown to
            # this version.
            query = query.filter(self.source.source.in_(list(self.sources.keys())))

        return query

    def _filter_by_date(self, query):
        """Filter query by observation date."""

        if self.start_date is not None:
            query = query.filter(Observation.mjd_stop >= self.start_date.mjd)

        if self.stop_date is not None:
            query = query.filter(Observation.mjd_start <= self.stop_date.mjd)

        return query

    def find_observations_containing_point(
        self, target: FixedTarget
    ) -> List[Observation]:
        """Find observations that contain the given point.


        Parameters
        ----------
        target : FixedTarget
            Position to search.


        Returns
        -------
        observations: list of Observation

        """

        terms: List[str] = self.indexer.query_point(target.ra.rad, target.dec.rad)

        q: Query = self.db.session.query(Observation)
        q = self._filter_by_source(q)
        q = self._filter_by_date(q)

        candidates: List[Observation] = q.filter(
            self.source.spatial_terms.overlap(terms)
        ).all()

        observations: List[Observation] = []
        for obs in candidates:
            fov_ra, fov_dec = core.polygon_string_to_arrays(obs.fov)
            if polygon_contains_point(fov_ra, fov_dec, target.ra.rad, target.dec.rad):
                observations.append(obs)
        if self.source != Observation:
            obsids: List[int] = [o.observation_id for o in observations]
            observations = (
                self.db.session.query(self.source).filter(
                    self.source.observation_id.in_(obsids)
                )
            ).all()

        return observations

    def find_observations_intersecting_cap(
        self, target: FixedTarget
    ) -> List[Observation]:
        """Find observations intersecting a spherical cap.

        The `padding` property is used for the cap radius, and
        `intersection_type` for the query.


        Parameters
        ----------
        target : FixedTarget
            The position to search.


        Returns
        -------
        observations : list of Observation

        """

        radius = np.radians(self.padding / 60)
        terms: List[str] = self.indexer.query_cap(target.ra.rad, target.dec.rad, radius)

        q: Query = self.db.session.query(Observation)
        q = self._filter_by_source(q)
        q = self._filter_by_date(q)

        q = q.filter(self.source.spatial_terms.overlap(terms))

        candidates: List[Observation] = q.all()

        observations: List[Observation] = []
        for obs in candidates:
            fov_ra, fov_dec = core.polygon_string_to_arrays(obs.fov)
            intersects = polygon_intersects_cap(
                fov_ra,
                fov_dec,
                target.ra.rad,
                target.dec.rad,
                radius,
                self.intersection_type.value,
            )
            if intersects:
                observations.append(obs)

        return observations

    def find_observations_intersecting_polygon(
        self, ra: np.ndarray, dec: np.ndarray
    ) -> List[Observation]:
        """Find observations intersecting given polygon.

        Parameters
        ----------
        ra, dec: array-like
            Polygon vertices in units of radians.  It is assumed that the
            area is < 2 pi steradians.

        Returns
        -------
        observations: list of Observation

        """

        _ra: np.ndarray = np.array(ra, float)
        _dec: np.ndarray = np.array(dec, float)
        terms: List[str] = self.indexer.query_polygon(_ra, _dec)

        q: Query = self.db.session.query(Observation)
        q = self._filter_by_source(q)
        q = self._filter_by_date(q)

        candidates: List[Observation] = q.filter(
            self.source.spatial_terms.overlap(terms)
        ).all()
        observations: List[Observation] = []
        for obs in candidates:
            fov_ra, fov_dec = core.polygon_string_to_arrays(obs.fov)
            if polygon_intersects_polygon(fov_ra, fov_dec, _ra, _dec):
                observations.append(obs)

        if self.source != Observation:
            obsids: List[int] = [o.observation_id for o in observations]
            observations = (
                self.db.session.query(self.source).filter(
                    self.source.observation_id.in_(obsids)
                )
            ).all()

        return observations

    def find_observations_intersecting_line(
        self,
        ra: np.ndarray,
        dec: np.ndarray,
        a: Optional[np.ndarray] = None,
        b: Optional[np.ndarray] = None,
        approximate: bool = False,
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
            _a = np.array(a, float)
            _b = np.array(b, float)
            if len(_b) != len(_a):
                raise ValueError("Size of a and b must match.")
            query_about = True

        observations: List[Observation] = []
        # search for the full line at once
        terms: List[str]
        if query_about:
            terms = self.indexer.query_about_line(_ra, _dec, _a, _b)[0]
        else:
            terms = self.indexer.query_line(_ra, _dec)

        q: Query = self.db.session.query(Observation)
        q = self._filter_by_source(q)
        q = self._filter_by_date(q)

        candidates: List[Observation] = q.filter(
            self.source.spatial_terms.overlap(terms)
        ).all()

        if not approximate:
            # test each observation for intersection with the observation FOV
            for obs in candidates:
                fov_ra, fov_dec = core.polygon_string_to_arrays(obs.fov)
                if query_about:
                    intersects = polygon_intersects_about_line(
                        fov_ra, fov_dec, _ra, _dec, _a, _b
                    )
                else:
                    intersects = polygon_intersects_line(fov_ra, fov_dec, _ra, _dec)
                if intersects:
                    observations.append(obs)
        else:
            observations = candidates

        if self.source != Observation:
            obsids: List[int] = [o.observation_id for o in observations]
            observations = (
                self.db.session.query(self.source).filter(
                    self.source.observation_id.in_(obsids)
                )
            ).all()

        return observations

    def find_observations_intersecting_line_at_time(
        self,
        ra: np.ndarray,
        dec: np.ndarray,
        mjd: np.ndarray,
        a: Optional[np.ndarray] = None,
        b: Optional[np.ndarray] = None,
        approximate: bool = False,
    ) -> List[Observation]:
        """Find observations intersecting given line at given times.

        Progress is published to ``self.search_logger``.


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
            raise ValueError("ra and dec must have same length")

        _a: Optional[np.ndarray] = None
        _b: Optional[np.ndarray] = None
        if a is not None or b is not None:
            _a: np.ndarray = np.array(a, float)
            _b: np.ndarray = np.array(b, float)
            if len(_a) != N or len(_b) != N:
                raise ValueError("ra, dec, a, and b must have same length")

        segments: Tuple[List[str], slice] = core.line_to_segment_query_terms(
            self.indexer,
            _ra,
            _dec,
            mjd,
            _a,
            _b,
            arc_limit=self.arc_limit,
            time_limit=self.time_limit,
        )

        self.logger.debug(
            "Splitting line of length %d into segments to test observation "
            "intersection at time.",
            N,
        )
        coords: SkyCoord = SkyCoord(_ra, _dec, unit="rad")
        arc_length: float = np.sum(coords[:-1].separation(coords[1:])).deg
        self.search_logger.info(
            "Searching coordinates with length %.1f deg and %.1f days",
            arc_length,
            np.ptp(mjd),
        )

        observations: List[Observation] = []
        n_segment_queries: int = 0
        n_matched_observations: int = 0

        terms: List[str]  # terms corresponding to each segment
        segment: slice  # slice corresponding to each segment
        for terms, segment in segments:
            n_segment_queries += 1
            q: Query = self.db.session.query(Observation)
            q = self._filter_by_source(q)

            # Only search based on mjd_start (alternatively, mjd_stop) so that
            # we can utilize a single database index.  But, add some padding to
            # the boundaries.  0.1 days = 2.4 hr... few survey integrations are
            # even that long.
            nearby_obs: List[Observation] = (
                q.filter(Observation.spatial_terms.overlap(terms))
                .filter(Observation.mjd_start >= (min(mjd[segment]) - 0.1))
                .filter(Observation.mjd_start <= max(mjd[segment]))
                .all()
            )

            # now do proper mjd_stop cut:
            nearby_obs = [
                obs for obs in nearby_obs if obs.mjd_stop >= min(mjd[segment])
            ]
            n_matched_observations += len(nearby_obs)

            # apply start_/stop_date range
            if self.start_date is not None:
                nearby_obs = [
                    obs for obs in nearby_obs if obs.mjd_start > self.start_date.mjd
                ]
            if self.stop_date is not None:
                nearby_obs = [
                    obs for obs in nearby_obs if obs.mjd_stop < self.stop_date.mjd
                ]

            if len(nearby_obs) > 0:
                if approximate:
                    observations.extend(nearby_obs)
                else:
                    # check for detailed intersection
                    self.logger.debug(
                        "Testing %d observations obtained between %s and %s "
                        "for detailed intersection.",
                        len(nearby_obs),
                        Time(min(mjd[segment]), format="mjd").iso[:10],
                        Time(max(mjd[segment]), format="mjd").iso[:10],
                    )
                    observations.extend(
                        core.test_line_intersection_with_observations_at_time(
                            nearby_obs,
                            _ra[segment],
                            _dec[segment],
                            mjd[segment],
                            None if a is None else _a[segment],
                            None if b is None else _b[segment],
                        )
                    )
        self.search_logger.info(
            "%d observation%s found",
            n_matched_observations if approximate else len(observations),
            "" if len(observations) == 1 else "s",
        )

        # duplicates can accumulate because each segment is searched
        # individually
        observations = list(set(observations))

        if self.source != Observation:
            obsids: List[int] = [o.observation_id for o in observations]
            observations = (
                self.db.session.query(self.source).filter(
                    self.source.observation_id.in_(obsids)
                )
            ).all()

        if approximate:
            self.logger.debug(
                "Tested %d segments, matched %d observations.",
                n_segment_queries,
                n_matched_observations,
            )
        else:
            self.logger.debug(
                "Tested %d segments, matched %d observations, %d intersections.",
                n_segment_queries,
                n_matched_observations,
                len(observations),
            )
        return observations

    def find_observations_by_ephemeris(
        self,
        eph: List[Ephemeris],
        approximate=False,
    ) -> List[Observation]:
        """Find observations covering given ephemeris.


        Parameters
        ----------
        eph: list of Ephemeris
            The ephemeris points to check.  Assumed to be continuous.

        approximate: bool, optional
            Do not check potential matches in detail.


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
            a, b = core.ephemeris_uncertainty_offsets(eph)
            obs = self.find_observations_intersecting_line_at_time(
                ra,
                dec,
                mjd,
                a=a + padding,
                b=b + padding,
                approximate=approximate,
            )
        elif self.padding > 0:
            obs = self.find_observations_intersecting_line_at_time(
                ra, dec, mjd, a=padding, b=padding, approximate=approximate
            )
        else:
            obs = self.find_observations_intersecting_line_at_time(
                ra, dec, mjd, approximate=approximate
            )

        self.logger.info(
            "%d observation%s found.", len(obs), "" if len(obs) == 1 else "s"
        )

        return obs
