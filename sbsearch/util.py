# Licensed with the 3-clause BSD license.  See LICENSE for details.
"""utility closet"""
import numpy as np
from astropy.time import Time
import astropy.coordinates as coords
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
    if len(constraints) > 0:
        expr, params = list(zip(*constraints))
        cmd = cmd + ' WHERE ' + ' AND '.join(expr)
        parameters.extend([p for p in params if p is not None])
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
    """Spherical interpolation by rotation.

    Parameters
    ----------
    c0, c1 : astropy.coordinates.SkyCoord
        RA and Dec coordinates of each point.

    t0, t1, t2 : float
        Time for each point (``t0``, ``t1``), and value to interpolate
        to (``t2``).


    Returns
    -------
    c2 : float
        Interpolated coordinate.

    """

    dt = (t2 - t0) / (t1 - t0)
    w = c0.separation(c1)

    a = sc2xyz(c0)
    b = sc2xyz(c1)
    n = np.cross(a, b)
    n /= np.sqrt((n**2).sum())

    c = vector_rotate(a, n, w * dt)
    d, dec, ra = coords.cartesian_to_spherical(*c)
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

    jda = jd[1] - half_step.to('day').value
    jdc = jd[1] + half_step.to('day').value
    a = spherical_interpolation(eph[0], eph[1], jd[0], jd[1], jda)
    b = eph[1]
    c = spherical_interpolation(eph[1], eph[2], jd[1], jd[2], jdc)
    x, y, z = list(zip(*[sc2xyz(sc) for sc in (a, b, c)]))
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


def vector_rotate(r, n, th):
    """Rotate vector `r` an angle `th` CCW about `n`.

    Parameters
    ----------
    r : array (3)
      Vector to rotate [x, y, z].
    n : array (3)
      Unit vector to rotate about.
    th : float or array
      The CCW angle to rotate by. [radians]

    Returns
    -------
    rp : ndarray
      The rotated vector [x, y, z].

    Notes
    -----
    Described in Goldstein p165, 2nd ed. Note that Goldstein presents
    the formula for clockwise rotation.

    """

    return (r * np.cos(-th) +
            n * (n * r).sum() * (1.0 - np.cos(-th)) +
            np.cross(r, n) * np.sin(-th))
