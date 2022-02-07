# Licensed with the 3-clause BSD license.  See LICENSE for details.
"""core

Helper functions for SBSearch.

"""

from typing import List, Tuple, Optional
import numpy as np

from sbsearch.exceptions import SBSException
from .model import Observation, Ephemeris
from .spatial import (  # pylint: disable=E0611
    SpatialIndexer,
    polygon_string_intersects_line,
    polygon_string_intersects_about_line,
)


def test_line_intersection_with_observations_at_time(
    obs: List[Observation],
    ra: np.ndarray,
    dec: np.ndarray,
    mjd: np.ndarray,
    a: Optional[np.ndarray] = None,
    b: Optional[np.ndarray] = None,
) -> List[Observation]:
    """Test each observation for intersection with the observation FOV.

    Comet and asteroid motion is non-linear, but this method uses a linear
    approximation.  (Non-linearity should be addressed with finer ephemeris
    steps.)  In order to minimize errors, only test the nearest segment(s) to
    the observation.  For example:

            0               1
    |----------|--------------------|

    Segment 0: t0 = 0 dt = 1 da = 10 deg

    Segment 1: t0 = 1 dt = 1 da = 20 deg

    Average proper motion: 30 deg / 2 days = 15 deg / day

    Linear interpolation to t = 1? --> 15 deg

    But we wanted 10 deg.

    """
    query_about: bool = a is not None
    _obs: List[Observation] = []
    if np.any(np.diff(mjd) <= 0):
        raise ValueError("Line segments must be monotonically increasing with time.")

    N: int = len(mjd)
    for o in obs:
        # find the nearest segment(s)
        i: int = max(np.searchsorted(mjd, o.mjd_start, side="right") - 1, 0)
        j: int = min(np.searchsorted(mjd, o.mjd_stop, side="right"), N - 1)
        segment: slice = slice(i, j + 1)

        dt = mjd[j] - mjd[i]
        line_start = (o.mjd_start - mjd[i]) / dt
        line_stop = (o.mjd_stop - mjd[i]) / dt
        try:
            if query_about:
                intersects = polygon_string_intersects_about_line(
                    o.fov,
                    ra[segment],
                    dec[segment],
                    a[segment],
                    b[segment],
                    line_start=line_start,
                    line_stop=line_stop,
                )
            else:
                intersects = polygon_string_intersects_line(
                    o.fov,
                    ra[segment],
                    dec[segment],
                    line_start=line_start,
                    line_stop=line_stop,
                )
        except ValueError as e:
            raise SBSException(
                f"Failed testing line intersection with {str(o)}."
            ) from e
        if intersects:
            _obs.append(o)

    return _obs


def line_to_segment_query_terms(
    indexer: SpatialIndexer,
    ra: np.ndarray,
    dec: np.ndarray,
    mjd: Optional[np.ndarray] = None,
    a: Optional[np.ndarray] = None,
    b: Optional[np.ndarray] = None,
    arc_limit: float = 0.17,
    time_limit: float = 365,
) -> Tuple[List[str], slice]:
    """Break-up a line, iterate over each segment's query terms.

    Query every ``arc_limit`` of motion (radians) or ``time_limit`` of
    observations (days), linear estimates are OK for this.  10 deg / 30 days =
    50 arcsec/hr.

    Inputs are not validated.

    ra, dec, a, b: ndarrays, radians.

    """

    N = len(ra)
    query_about: bool = a is not None
    time_test: bool = mjd is not None
    da: float = 0
    i: int = 0  # first index to test
    j: int = 1  # last
    while True:
        segment: slice = slice(i, j + 1)
        da += np.hypot(ra[j - 1] - ra[j], dec[j - 1] - dec[j])
        dt: float = mjd[j] - mjd[i] if time_test else 0

        if (da < arc_limit) and (dt < time_limit) and (j < N - 1):
            j += 1
            continue

        terms: List[str]
        if query_about:
            terms = indexer.query_about_line(
                ra[segment], dec[segment], a[segment], b[segment]
            )[0]
        else:
            terms = indexer.query_line(ra[segment], dec[segment])

        yield terms, segment

        da = 0
        i = j
        j += 1

        if i >= N - 1:
            break


def ephemeris_uncertainty_offsets(eph: List[Ephemeris]) -> Tuple[np.ndarray]:
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
        raise ValueError("Must have at least 2 ephemeris points.")

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
    P: np.ndarray = np.array((np.cos(pa + np.pi / 2), np.sin(pa + np.pi / 2)))

    # setup uncertainty vectors
    A: np.ndarray = unc_a * np.array((np.cos(unc_theta), np.sin(unc_theta)))
    B: np.ndarray = unc_b * np.array(
        (np.cos(unc_theta + np.pi / 2), np.sin(unc_theta + np.pi / 2))
    )

    # compute offsets a, b
    a = np.max(np.abs((np.sum(A * R, 0), np.sum(B * R, 0))), 0)
    b = np.max(np.abs((np.sum(A * P, 0), np.sum(B * P, 0))), 0)

    return a, b
