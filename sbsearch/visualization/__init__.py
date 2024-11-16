# Licensed with the 3-clause BSD license.  See LICENSE for details.

__all__ = [
    "vertices_to_spherical_polygon",
    "observation_to_spherical_polygon",
    "term_to_spherical_polygon",
    "plot_terms",
    "plot_observations",
    "plot_fixed_target",
]

from typing import Iterable, List
import numpy as np

try:
    import matplotlib as mpl
    from matplotlib.patches import Polygon, Rectangle
    from matplotlib.collections import PatchCollection
    from spherical_geometry.polygon import (
        SingleSphericalPolygon,
        great_circle_arc,
        vector,
    )
except ImportError:  # pragma: no cover
    mpl = None

import astropy.units as u
from astropy.visualization.wcsaxes import WCSAxes


from ..core import polygon_string_to_arrays
from ..model import Observation
from ..spatial import term_to_cell_vertices


def vertices_to_spherical_polygon(ra: np.ndarray, dec: np.ndarray, **kwargs) -> Polygon:
    """Convert a set of vertices a spherical polygon.

    The edges of a spherical polygon are great circles.


    Parameters
    ----------
    ra, dec : array
        Right ascension and declination as in units of radians.

    kwargs :
        Any `matplotlib.patches.Polygon` keyword argument.


    Returns
    -------
    poly : `matplotlib.patches.Polygon`

    """

    if mpl is None:  # pragma: no cover
        raise RuntimeError("Requires matplotlib and spherical_geometry")

    polygon: SingleSphericalPolygon = SingleSphericalPolygon.from_radec(
        ra, dec, degrees=False
    )

    _ra: List[float]
    _dec: List[float]
    _ra, _dec = [], []
    for A, B in zip(polygon.points[0:-1], polygon.points[1:]):
        length: float = great_circle_arc.length(A, B)
        if not np.isfinite(length):
            length = 2

        interpolated: np.ndarray = great_circle_arc.interpolate(A, B, length * 4)

        lon: np.ndarray
        lat: np.ndarray
        lon, lat = vector.vector_to_lonlat(
            interpolated[:, 0],
            interpolated[:, 1],
            interpolated[:, 2],
            degrees=False,
        )

        lon0: float
        lat0: float
        lon1: float
        lat1: float
        for lon0, lat0, lon1, lat1 in zip(lon[0:-1], lat[0:-1], lon[1:], lat[1:]):
            _ra.append(lon0)
            _dec.append(lat0)

    if len(_ra) > 0:
        _ra.append(lon1)
        _dec.append(lat1)

    vertices: u.Quantity = np.degrees([_ra, _dec]).transpose() * u.deg
    return Polygon(vertices, **kwargs)


def observation_to_spherical_polygon(observation: Observation, **kwargs) -> Polygon:
    """Convert an observation object to a spherical polygon.

    The edges of a spherical polygon are great circles.


    Parameters
    ----------
    observation : list of Observation
        The fov property will be converted.

    kwargs :
        Any `matplotlib.patches.Polygon` keyword argument.


    Returns
    -------
    poly : `matplotlib.patches.Polygon`

    """

    return vertices_to_spherical_polygon(
        *polygon_string_to_arrays(observation.fov), **kwargs
    )


def term_to_spherical_polygon(term: str, **kwargs) -> Polygon:
    """Convert S2 cell to spherical polygon.

    The edges of a spherical polygon are great circles.


    Parameters
    ----------
    term : str
        S2 query/index cell term/token.

    kwargs :
        Any `matplotlib.patches.Polygon` keyword argument.


    Returns
    -------
    poly : `matplotlib.patches.Polygon`

    """

    ra: np.ndarray
    dec: np.ndarray
    ra, dec = term_to_cell_vertices(term.lstrip("$"))
    ra[ra < 0] += 2 * np.pi

    return vertices_to_spherical_polygon(ra, dec, **kwargs)


def plot_polygons(ax: WCSAxes, polygons: Iterable, **kwargs) -> Rectangle:
    """Plot polygons as a `PatchCollection`.


    Parameters
    ----------

    ax : `astropy.visualization.wcsaxes.WCSAxes`
        The axes to plot to.

    polygons : array
        The N polygons with M vertices to plot as an (N, M, 2) array.  The
        coordinates are RA and Dec in degrees.

    **kwargs
        Keyword arguments are passed to the `PatchCollection` object.


    Returns
    -------
    handle : Rectangle
        A rectangle with the same style to use as a legend handle.

    """

    collection = PatchCollection(
        polygons,
        transform=ax.get_transform("icrs"),
        **kwargs,
    )

    ax.add_artist(collection)
    return Rectangle((0, 0), 1, 1, **kwargs)


def plot_terms(ax: WCSAxes, terms: List[str], **kwargs) -> Rectangle:
    """Plot S2 cell terms/tokens.


    Parameters
    ----------

    ax : `astropy.visualization.wcsaxes.WCSAxes`
        The axes to plot to.

    terms : list of str
        The terms to plot.

    **kwargs
        Keyword arguments are passed to the `PatchCollection` object.


    Returns
    -------
    handle : Rectangle
        A rectangle with the same style to use as a legend handle.

    """

    return plot_polygons(
        ax,
        [term_to_spherical_polygon(term) for term in terms],
        **kwargs,
    )


def plot_observations(
    ax: WCSAxes, observations: List[Observation], **kwargs
) -> Rectangle:
    """Plot observation fields of view.


    Parameters
    ----------

    ax : `astropy.visualization.wcsaxes.WCSAxes`
        The axes to plot to.

    observations : list of Observation
        The observations to plot.

    **kwargs
        Keyword arguments are passed to the `PatchCollection` object.


    Returns
    -------
    handle : Rectangle
        A rectangle with the same style to use as a legend handle.

    """

    return plot_polygons(
        ax,
        [observation_to_spherical_polygon(obs) for obs in observations],
        **kwargs,
    )
