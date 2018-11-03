# Licensed with the 3-clause BSD license.  See LICENSE for details.
"""utility closet"""
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord


def assemble_sql(cmd, parameters, constraints):
    """Assemble a SQL statement.

    Parameters
    ----------
    cmd : string
        Left-hand side of the statement.

    parameters : list
        Parameters for substitution (via sqlite3 parameter
        substitution).

    constraints : list of tuple
        Each constraint is a SQL expression and an optional parameter
        for subsitution into the expression.  If no parameter is used,
        set it to ``None``.

    """
    if len(contraints) > 0:
        expr, params = list(zip(constraints))
        cmd = cmd + ' WHERE ' + ' AND '.join(expr)
        parameters.extend(params)
    return cmd, parameters


def date_constraints(jd_start, jd_stop):
    """Add date constraints for assemble_sql()."""
    constraints = []
    if jd_start is not None:
        constraints.append(('jd>=?', jd_start))

    if jd_stop is not None:
        constraints.append(('jd<=?', jd_stop))

    return constraints


def iterate_over(cursor):
    while True:
        rows = cursor.fetchmany()
        if not rows:
            return
        for row in rows:
            yield row


def interior_test(self, point, corners):
    """Test if point is interior to corners assuming spherical geometry."""
    pass


def rd2xyz(ra, dec):
    """RA, Dec (radians or Angle) to Cartesian coordinates."""
    return np.array((np.cos(dec) * np.cos(ra),
                     np.cos(dec) * np.sin(ra),
                     np.sin(dec)))


def sc2xyz(sc):
    """SkyCoord to Cartesian coordinates."""
    return np.array((np.cos(sc.dec) * np.cos(sc.ra),
                     np.cos(sc.dec) * np.sin(sc.ra),
                     np.sin(sc.dec)))


def spherical_interpolation(c0, c1, t0, t1, t2):
    """Spherical interpolation between two points.

    Parameters
    ----------
    c0, c1 : astropy.coordinates.SkyCoord
        RA and Dec coordinates of each point.

    t0, t1, t2 : float
        Time for each point (``t1``, ``t2``), and the result
        (``t2``).


    Returns
    -------
    c2 : float
        Interpolated coordinate.

    """

    dt = (t2 - t0) / (t1 - t0)
    w = c0.separation(c1)

    p1 = np.sin((1 - dt) * w.rad) / np.sin(w.rad)
    p2 = np.sin(dt * w.rad) / np.sin(w.rad)

    ra = p1 * c1.ra + p2 * c2.ra
    dec = p1 * c1.dec + p2 * c2.dec

    return SkyCoord(ra, dec)


def eph_to_limits(jd, eph, half_step):
    """Specialized for the ephemeris R-tree.

    Take a 3-point ephemeris and find the x, y, z, and t range that is
    centered on the second point, with a length of ``half_step * 2``.

    Parameters
    ----------
    jd : array
        Julian-date of points to interpolate between.

    eph : SkyCoord
        RA, Dec.

    half_step : astropy.units.Quantity
        Half the step size between points in days.

    """

    jda = jd - half_step
    jdc = jd + half_step
    a = util.spherical_interpolation(
        eph[0], eph[1], jd[0], jd[1], jd - half_step)
    b = eph[1]
    c = util.spherical_interpolation(
        eph[1], eph[2], jd[1], jd[2], jd + half_step)
    x, y, z = list(zip([sc2xyz(sc) for sc in (a, b, c)]))
    return jda, jdc, min(x), max(x), min(y), max(y), min(z), max(z)


def epochs_to_time(epochs, scale='utc'):
    """Flexible time input to `~astropy.time.Time` object.

    Parameters
    ----------
    epochs : iteratable
        May be integers or floats for Julian date, or any object
        parseable by `~astropy.time.Time`.

    scale : string, optional
        Time scale.

    Returns
    -------
    times : `~astropy.time.Time`

    """

    times = []
    for epoch in epochs:
        if isinstance(epoch, (float, int)):
            format = 'jd'
        else:
            format = None

        times.append(Time(epoch, format=format, scale=scale))

    return Time(times)
