# Licensed with the 3-clause BSD license.  See LICENSE for details.
"""utility closet"""
import struct
import numpy as np
from astropy.time import Time
import astropy.coordinates as coords
from astropy.coordinates import Angle
from astropy.coordinates.angle_utilities import angular_separation
import astropy.units as u


class RADec:
    def __init__(self, *args, unit=None):
        if isinstance(args[0], RADec):
            self.ra = args[0].ra
            self.dec = args[0].dec
        else:
            if np.iterable(args[0]) and len(args) == 1:
                a = np.array(args[0])
                ra = a[..., 0]
                dec = a[..., 1]
            elif len(args) == 2:
                ra = args[0]
                dec = args[1]
            else:
                raise ValueError('Unknown input: {}'.format(args))

            self.ra = Angle(ra, unit=unit).rad
            self.dec = Angle(dec, unit=unit).rad

    def __repr__(self):
        return "<RADec: ra={}, dec={}>".format(self.ra, self.dec)

    def __len__(self):
        return np.size(self.ra)

    def __getitem__(self, k):
        return RADec(self.ra[k], self.dec[k], unit='rad')

    def separation(self, other):
        return angular_separation(self.ra, self.dec, other.ra, other.dec)

    @property
    def xyz(self):
        return rd2xyz(self.ra, self.dec)


def assemble_sql(cmd, parameters, constraints, inner_join=None):
    """Assemble a SQL statement.

    Parameters
    ----------
    cmd : string
        Left-hand side of the statement.

    parameters : list
        Parameters for substitution (via sqlite3 parameter
        substitution).

    constraints : list of tuple
        Each constraint is a SQL expression and optional parameters
        for subsitution into the expression.  If no parameter is used,
        set it to ``None``.

    inner_join : string or list of strings, optional
        List of inner join constraints, e.g., `'obs USING (obsid)'`.

    """
    if isinstance(inner_join, str):
        inner_join = [inner_join]

    if inner_join:
        for ij in inner_join:
            cmd += ' INNER JOIN ' + ij

    if len(constraints) > 0:
        expr, params = list(zip(*constraints))
        cmd = cmd + ' WHERE ' + ' AND '.join(expr)
        for p in params:
            if p is None:
                continue
            if isinstance(p, (list, tuple)):
                parameters.extend(p)
            else:
                parameters.append(p)
    return cmd, parameters


def date_constraints(start, stop, column='jd'):
    """Add date constraints for assemble_sql()."""
    constraints = []
    if start is None:
        jd_start = None
    else:
        if isinstance(start, (float, int)):
            jd_start = start
        else:
            jd_start = Time(start).jd

    if stop is None:
        jd_stop = None
    else:
        if isinstance(stop, (float, int)):
            jd_stop = stop
        else:
            jd_stop = Time(stop).jd

    if jd_start is not None:
        constraints.append((column + '>=?', jd_start))

    if jd_stop is not None:
        constraints.append((column + '<=?', jd_stop))

    return constraints


def eph_to_limits(eph, jd, half_step):
    """Specialized for the ephemeris R-tree.

    Take a 3-point ephemeris and find the x, y, z, and t range that is
    centered on the second point, with a length of ``half_step * 2``.

    Parameters
    ----------
    eph : RADec
        RA, Dec.

    jd : array
        Julian-date of points to interpolate between.

    half_step : astropy.units.Quantity
        Half the step size between points in days.

    """

    dt = u.Quantity(half_step, 'day').value
    mjd = np.array(jd) - 2400000.5
    mjda = mjd[1] - dt
    mjdc = mjd[1] + dt

    if np.allclose((eph[0].ra, eph[0].dec), (eph[1].ra, eph[1].dec)):
        a = eph[0]
    else:
        a = spherical_interpolation(eph[0], eph[1], mjd[0], mjd[1], mjda)

    b = eph[1]

    if np.allclose((eph[1].ra, eph[1].dec), (eph[2].ra, eph[2].dec)):
        c = eph[2]
    else:
        c = spherical_interpolation(eph[1], eph[2], mjd[1], mjd[2], mjdc)

    x, y, z = list(zip(*[sc.xyz for sc in (a, b, c)]))
    return {
        'mjd0': mjda,
        'mjd1': mjdc,
        'x0': min(x),
        'x1': max(x),
        'y0': min(y),
        'y1': max(y),
        'z0': min(z),
        'z1': max(z)
    }


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


def epochs_to_jd(epochs):
    """Flexible time input to Julian date.

    Parameters
    ----------
    epochs : iteratable
        May be integers or floats for Julian date, or any object
        parseable by `~astropy.time.Time`.  ``None`` items are left
        as-is.

    Returns
    -------
    jd : list

    """

    jd = []
    for epoch in epochs:
        if isinstance(epoch, (float, int)) or epoch is None:
            jd.append(epoch)
        else:
            jd.append(Time(epoch).jd)

    return jd


def fov2points(fov):
    """Convert from obs database FOV to arrays of points.

    Parameters
    ----------
    fov : bytes array
        FOV blob from the database.

    Returns
    -------
    ra, dec : ndarray
        RA and Dec arrays (radians).

    """
    N = len(fov) // 80
    x = struct.unpack('{}d'.format(10 * N), fov)
    ra = np.array(x[::2])
    dec = np.array(x[1::2])
    return ra, dec


def iterate_over(cursor):
    """Iterate over SQLite cursour via fetchmany."""
    while True:
        rows = cursor.fetchmany()
        if not rows:
            return
        for row in rows:
            yield row


def rd2xyz(ra, dec):
    """RA, Dec (radians or Angle) to Cartesian coordinates."""
    return np.array((np.cos(dec) * np.cos(ra),
                     np.cos(dec) * np.sin(ra),
                     np.sin(dec)))


def spherical_interpolation(c0, c1, t0, t1, t2):
    """Spherical interpolation by rotation.


    Parameters
    ----------
    c0, c1 : RADec
        Coordinates of each point.

    t0, t1, t2 : float
        Time for each point (``t0``, ``t1``), and value to interpolate
        to (``t2``).


    Returns
    -------
    c2 : RADec
        Interpolated coordinate.

    """

    if t0 == t1:
        raise ValueError('t0 == t1')

    if t2 == t0:
        return c0

    if t2 == t1:
        return c1

    dt = (t2 - t0) / (t1 - t0)
    w = c0.separation(c1)

    a = c0.xyz
    b = c1.xyz
    n = np.cross(a, b)
    n /= np.sqrt((n**2).sum())

    c = vector_rotate(a, n, w * dt)
    d, dec, ra = coords.cartesian_to_spherical(*c)
    return RADec(ra, dec)


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


def vmag_from_eph(eph, ignore_zero=True, missing=99):
    """Get most relevant magnitude estimate from ephemeris.

    Parameters
    ----------
    eph : `~sbpy.data.Ephem`
        Ephemeris.

    ignore_zero : bool, optional
        ``Ephem`` does not support masking, so ephemerides are
        populated with zeros.  Set to ``True`` to ignore them and use
        another magnitude estimate, if available.  Order of
        preference: Tmag, Nmag, V

    missing : float, optional
        Use this value for missing magnitudes.

    """

    vmag = missing * np.ones(len(eph.table))
    if ignore_zero:
        for k in ['V', 'Nmag', 'Tmag']:
            if k in eph.table.colnames:
                i = eph[k].value != 0
                vmag[i] = eph[k][i].value
    else:
        for k in ['Tmag', 'Nmag', 'V']:
            if k in eph.table.colnames:
                vmag = eph[k].value
                break
    return vmag
