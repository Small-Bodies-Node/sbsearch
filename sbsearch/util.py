# Licensed with the 3-clause BSD license.  See LICENSE for details.
"""utility closet"""
import struct
import numpy as np
from astropy.time import Time
import astropy.coordinates as coords
from astropy.coordinates import Angle
from astropy.coordinates.angle_utilities import angular_separation
import astropy.units as u
import geoalchemy2
from . import schema


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

            self.ra = Angle(ra, unit=unit)
            self.dec = Angle(dec, unit=unit)

    @classmethod
    def from_eph(self, eph):
        """Initialize from Eph or list of Eph."""

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


class FieldOfView:
    """Polygon on the sky.

    Parameters
    ----------
    vertices : RADec
        FOV corners, in order.

    """

    def __init__(self, vertices):
        if not isinstance(vertices, RADec):
            raise TypeError('vertices must be RADec')
        self.vertices = vertices

    def __str__(self):
        """PostGIS formatted string."""
        vertices = [v for v in self.vertices] + [self.vertices[0]]
        polygon = ','.join(['{} {}'.format(float(v.ra.deg), float(v.dec.deg))
                            for v in vertices])
        return 'SRID=40001;POLYGON(({}))'.format(polygon)


class Line:
    """Line on the sky.

    Parameters
    ----------
    vertices : RADec
        Line vertices, must have at least 2 points.

    """

    def __init__(self, vertices):
        self.vertices = vertices

    @classmethod
    def from_eph(cls, eph):
        """Initialize line from Eph object.

        Returns
        -------
        line : Line
            Line representing ``eph.segment``.

        """
        g = geoalchemy2.functions.ST_AsGeoJSON(eph.segment)
        coords = np.array(g['geometry']['coordinates'])
        return super()(RADec(coords, unit='deg'))

    def __str__(self):
        """PostGIS formatted string."""
        vertices = [v for v in self.vertices]
        polygon = ','.join(['{} {}'.format(v.ra.deg, v.dec.deg)
                            for v in vertices])
        return 'SRID=40001;LINESTRING(({}))'.format(vertices)


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


def filter_by_date_range(query, start, stop, column):
    """Filter SQLAlchemy query by date range.


    Parameters
    ----------

    query : sqlalchemy Query
        The query to filter.

    start, stop : int, float, str, None
        Integer or float for Julian date, else a UTC string parseable
        by `~astropy.time.Time`.  Use ``None`` for no limit.

    column : sqlalchemy Column
        Filter this column.


    Returns
    -------
    revised_query

    """

    if start is not None:
        if isinstance(start, str):
            start = Time(start).jd

        query = query.filter(column >= start)

    if stop is not None:
        if isinstance(stop, str):
            stop = Time(start).jd

        query = query.filter(column <= stop)

    return query


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
