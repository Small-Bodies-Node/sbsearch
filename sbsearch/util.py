# Licensed with the 3-clause BSD license.  See LICENSE for details.
"""utility closet"""
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
                ra = a[:, 0]
                dec = a[:, 1]
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

    if eph[0] == eph[1]:
        a = eph[0]
    else:
        a = spherical_interpolation(eph[0], eph[1], mjd[0], mjd[1], mjda)

    b = eph[1]

    if eph[1] == eph[2]:
        c = eph[2]
    else:
        c = spherical_interpolation(eph[1], eph[2], mjd[1], mjd[2], mjdc)

    x, y, z = list(zip(*[sc.xyz for sc in (a, b, c)]))
    return mjda, mjdc, min(x), max(x), min(y), max(y), min(z), max(z)


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


def interior_test(point, corners):
    """Test if point is interior to corners assuming spherical geometry.


    Parameters
    ----------
    point : RADec
        Point to test.

    corners : RADec
        Points describing a spherical rectangle.


    Returns
    -------
    interior : 
        ``True`` if the point falls inside the rectangle.

    """

    # 0, k and i, j are opposite corners
    i, j, k = corners[0].separation(corners[1:]).argsort() + 1
    segments = np.array((
        corners[0].separation(corners[i]),
        corners[i].separation(corners[k]),
        corners[k].separation(corners[j]),
        corners[j].separation(corners[0]),
        corners[j].separation(corners[i])
    ))

    # measure area of the rectangle
    area_r = spherical_triangle_area(segments[0], segments[1], segments[4])
    area_r += spherical_triangle_area(segments[2], segments[3], segments[4])

    # measure area of all triangles made by point and rectangle segments
    d = point.separation(corners)
    area_p = spherical_triangle_area(segments[0], d[0], d[i])
    area_p += spherical_triangle_area(segments[1], d[i], d[k])
    area_p += spherical_triangle_area(segments[2], d[k], d[j])
    area_p += spherical_triangle_area(segments[3], d[j], d[0])

    return np.isclose(area_r, area_p, rtol=1e-5)


def interior_test1(point, corners):
    """Test if point is interior to corners assuming spherical geometry.


    Parameters
    ----------
    point : `~astropy.coordinates.SkyCoord`
        Point to test.

    corners : `~astropy.coordinates.SkyCoord`
        Points describing a spherical rectangle.


    Returns
    -------
    interior : 
        ``True`` if the point falls inside the rectangle.

    """

    pi = np.pi

    # 0, k and i, j are opposite corners
    i, j, k = corners[0].separation(corners[1:]).argsort() + 1
    segments = np.array((
        corners[0].separation(corners[i]).rad,
        corners[i].separation(corners[k]).rad,
        corners[k].separation(corners[j]).rad,
        corners[j].separation(corners[0]).rad,
        corners[j].separation(corners[i]).rad
    ))

    # measure area of the rectangle
    area_r = spherical_triangle_area(segments[0], segments[1], segments[4])
    area_r += spherical_triangle_area(segments[2], segments[3], segments[4])

    # measure area of all triangles made by point and rectangle segments
    d = point.separation(corners).rad
    area_p = spherical_triangle_area(segments[0], d[0], d[i])
    area_p += spherical_triangle_area(segments[1], d[i], d[k])
    area_p += spherical_triangle_area(segments[2], d[k], d[j])
    area_p += spherical_triangle_area(segments[3], d[j], d[0])

    return np.isclose(area_r, area_p, rtol=1e-5)


def interior_test2(point, corners):
    """Test if point is interior to corners assuming spherical geometry.


    Parameters
    ----------
    point : `~astropy.coordinates.SkyCoord`
        Point to test.

    corners : `~astropy.coordinates.SkyCoord`
        Points describing a spherical rectangle.


    Returns
    -------
    interior : 
        ``True`` if the point falls inside the rectangle.

    """

    # first three segments make 1/2 the rectangle, next four are
    # point to each rectangle corner
    segments = np.empty(9)

    # 0, k and i, j are opposite corners
    i, j, k = corners[0].separation(corners[1:]).argsort() + 1
    segments[:2] = corners[i].separation(corners[[0, k]]).rad
    segments[2:4] = corners[k].separation(corners[[i, j]]).rad
    segments[4] = corners[j].separation(corners[0]).rad
    segments[5:] = point.separation(corners[[0, i, j, k]])

    # make the triangles
    tri = np.empty((6, 3))
    tri[0] = segments[[0, 2, 1]]
    tri[1] = segments[[1, 3, 4]]
    tri[2] = segments[[0, 5, 6]]
    tri[3] = segments[[2, 6, 8]]
    tri[4] = segments[[3, 7, 8]]
    tri[5] = segments[[4, 5, 7]]

    # measure area of the rectangle
    ca, cb, cc = np.cos((tri[:, 0], tri[:, 1], tri[:, 2]))
    sa, sb, sc = np.sin((tri[:, 0], tri[:, 1], tri[:, 2]))
    A = np.arccos((ca - cb * cc) / (sb * sc))
    B = np.arccos((cb - cc * ca) / (sc * sa))
    C = np.arccos((cc - ca * cb) / (sa * sb))

    area = A + B + C - np.pi
    area_r = area[:2].sum()
    area_p = area[2:].sum()

    return np.isclose(area_r, area_p, rtol=1e-5)


def interior_test3(point, corners):
    """Test if point is interior to corners assuming spherical geometry.


    Parameters
    ----------
    point : `~astropy.coordinates.SkyCoord`
        Point to test.

    corners : `~astropy.coordinates.SkyCoord`
        Points describing a spherical rectangle.


    Returns
    -------
    interior : 
        ``True`` if the point falls inside the rectangle.

    """

    # first three segments make 1/2 the rectangle, next four are
    # point to each rectangle corner
    segments = np.empty(9)

    # 0, k and i, j are opposite corners
    i, j, k = corners[0].separation(corners[1:]).argsort() + 1
    segments[:2] = corners[i].separation(corners[[0, k]]).rad
    segments[2:4] = corners[k].separation(corners[[i, j]]).rad
    segments[4] = corners[j].separation(corners[0]).rad
    segments[5:] = point.separation(corners[[0, i, j, k]])

    # make the triangles
    tri = np.empty((6, 3))
    tri[0] = segments[[0, 2, 1]]
    tri[1] = segments[[1, 3, 4]]
    tri[2] = segments[[0, 5, 6]]
    tri[3] = segments[[2, 6, 8]]
    tri[4] = segments[[3, 7, 8]]
    tri[5] = segments[[4, 5, 7]]

    # measure area of the rectangle
    ca, cb, cc = np.cos((tri[:, 0], tri[:, 1], tri[:, 2]))
    sa, sb = np.sin((tri[:, 0], tri[:, 1]))
    ta2, tb2 = np.tan((tri[:, 0] / 2, tri[:, 1] / 2))
    cosC = (cc - ca * cb) / (sa * sb)
    sinC = np.sqrt(1 - cosC**2)

    tab2 = ta2 * tb2
    area = 2 * np.arctan(tab2 * sinC / (1 + tab2 * cosC))
    area_r = area[:2].sum()
    area_p = area[2:].sum()

    return np.isclose(area_r, area_p, rtol=1e-5)


def interior_test4(point, corners):
    """Test if point is interior to corners assuming spherical geometry.


    Parameters
    ----------
    point : array
        RA, Dec of point to test, radians.

    corners : array
        Corners (RA, Dec) of rectangle to test, 4x2, radians.


    Returns
    -------
    interior : 
        ``True`` if the point falls inside the rectangle.

    """

    # first three segments make 1/2 the rectangle, next four are
    # point to each rectangle corner
    segments = np.empty(9)

    corners = np.array(corners)

    # 0, k and i, j are opposite corners
    d = angular_separation(corners[0, 0], corners[0, 1],
                           corners[1:, 0], corners[1:, 1])
    i, j, k = np.argsort(d) + 1
    segments[:2] = angular_separation(
        corners[i, 0], corners[i, 1], corners[[0, k], 0], corners[[0, k], 1])
    segments[2:4] = angular_separation(
        corners[k, 0], corners[k, 1], corners[[i, j], 0], corners[[i, j], 1])
    segments[4] = d[j - 1]
    segments[5:] = angular_separation(
        point[0], point[1], corners[:, 0], corners[:, 1])

    # make the triangles
    tri = np.empty((6, 3))
    tri[0] = segments[[0, 2, 1]]
    tri[1] = segments[[1, 3, 4]]
    tri[2] = segments[[0, 5, 6]]
    tri[3] = segments[[2, 6, 8]]
    tri[4] = segments[[3, 7, 8]]
    tri[5] = segments[[4, 5, 7]]

    # measure area of the rectangle
    ca, cb, cc = np.cos((tri[:, 0], tri[:, 1], tri[:, 2]))
    sa, sb = np.sin((tri[:, 0], tri[:, 1]))
    ta2, tb2 = np.tan((tri[:, 0] / 2, tri[:, 1] / 2))
    cosC = (cc - ca * cb) / (sa * sb)
    sinC = np.sqrt(1 - cosC**2)

    tab2 = ta2 * tb2
    area = 2 * np.arctan(tab2 * sinC / (1 + tab2 * cosC))
    area_r = area[:2].sum()
    area_p = area[2:].sum()

    return np.isclose(area_r, area_p, rtol=1e-5)


def interior_test5(point, corners):
    """Test if point is interior to corners assuming spherical geometry.


    Parameters
    ----------
    point : `~astropy.coordinates.SkyCoord`
        Point to test.

    corners : `~astropy.coordinates.SkyCoord`
        Points describing a spherical rectangle.


    Returns
    -------
    interior : 
        ``True`` if the point falls inside the rectangle.

    """

    #p = point.ra.rad, point.dec.rad
    #r = np.c_[corners.ra.rad, corners.dec.rad]
    # return interior_test4(p, r)

    # first three segments make 1/2 the rectangle, next four are
    # point to each rectangle corner
    segments = np.empty(9)

    _point = (point.ra.rad, point.dec.rad)
    _corners = np.c_[corners.ra.rad, corners.dec.rad]

    # 0, k and i, j are opposite corners
    d = angular_separation(_corners[0, 0], _corners[0, 1],
                           _corners[1:, 0], _corners[1:, 1])
    i, j, k = np.argsort(d) + 1
    segments[:2] = angular_separation(
        _corners[i, 0], _corners[i, 1],
        _corners[[0, k], 0], _corners[[0, k], 1])
    segments[2:4] = angular_separation(
        _corners[k, 0], _corners[k, 1],
        _corners[[i, j], 0], _corners[[i, j], 1])
    segments[4] = d[j - 1]
    segments[5:] = angular_separation(
        _point[0], _point[1], _corners[:, 0], _corners[:, 1])

    # make the triangles
    tri = np.empty((6, 3))
    tri[0] = segments[[0, 2, 1]]
    tri[1] = segments[[1, 3, 4]]
    tri[2] = segments[[0, 5, 6]]
    tri[3] = segments[[2, 6, 8]]
    tri[4] = segments[[3, 7, 8]]
    tri[5] = segments[[4, 5, 7]]

    # measure area of the rectangle
    ca, cb, cc = np.cos((tri[:, 0], tri[:, 1], tri[:, 2]))
    sa, sb = np.sin((tri[:, 0], tri[:, 1]))
    ta2, tb2 = np.tan((tri[:, 0] / 2, tri[:, 1] / 2))
    cosC = (cc - ca * cb) / (sa * sb)
    sinC = np.sqrt(1 - cosC**2)

    tab2 = ta2 * tb2
    area = 2 * np.arctan(tab2 * sinC / (1 + tab2 * cosC))
    area_r = area[:2].sum()
    area_p = area[2:].sum()

    return np.isclose(area_r, area_p, rtol=1e-5)


def interior_test6(point, corners):
    """Test if point is interior to corners assuming spherical geometry.


    Parameters
    ----------
    point : RADec
        RA, Dec of point to test.

    corners : RADec
        Corners of rectangle to test.


    Returns
    -------
    interior : 
        ``True`` if the point falls inside the rectangle.

    """

    # first three segments make 1/2 the rectangle, next four are
    # point to each rectangle corner
    segments = np.empty(9)

    # 0, k and i, j are opposite corners
    d = corners[0].separation(corners[1:])
    i, j, k = np.argsort(d) + 1
    segments[:2] = corners[i].separation(corners[[0, k]])
    segments[2:4] = corners[k].separation(corners[[i, j]])
    segments[4] = d[j - 1]
    segments[5:] = point.separation(corners[[0, i, j, k]])

    # make the triangles
    tri = np.empty((6, 3))
    tri[0] = segments[[0, 2, 1]]
    tri[1] = segments[[1, 3, 4]]
    tri[2] = segments[[0, 5, 6]]
    tri[3] = segments[[2, 6, 8]]
    tri[4] = segments[[3, 7, 8]]
    tri[5] = segments[[4, 5, 7]]

    # measure area of the rectangle
    ca, cb, cc = np.cos((tri[:, 0], tri[:, 1], tri[:, 2]))
    sa, sb = np.sin((tri[:, 0], tri[:, 1]))
    ta2, tb2 = np.tan((tri[:, 0] / 2, tri[:, 1] / 2))
    cosC = (cc - ca * cb) / (sa * sb)
    sinC = np.sqrt(1 - cosC**2)

    tab2 = ta2 * tb2
    area = 2 * np.arctan(tab2 * sinC / (1 + tab2 * cosC))
    area_r = area[:2].sum()
    area_p = area[2:].sum()

    return np.isclose(area_r, area_p, rtol=1e-5)


def spherical_triangle_area(a, b, c):
    """Area of spherical triangle."""
    ca, cb, cc = np.cos((a, b, c))
    sa, sb, sc = np.sin((a, b, c))
    A = np.arccos((ca - cb * cc) / (sb * sc))
    B = np.arccos((cb - cc * ca) / (sc * sa))
    C = np.arccos((cc - ca * cb) / (sa * sb))

    return A + B + C - np.pi


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
        return c0

    if t2 == t0:
        return c0

    if t1 == t0:
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

    vmag = missing * np.ones(len(eph))
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
