"""
plot S2 text formatted polygons
"""

import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import TABLEAU_COLORS
import astropy.units as u
from astropy.coordinates import offset_by
import s2geometry as s2
from mskpy import minmax

# plot styles
styles = {
    "query point": dict(
        zorder=4, marker="*", color="y", edgecolor="k", linewidths=0.5, s=50
    ),
    "query area": dict(
        zorder=3, lw=1, color="tab:blue", fc=TABLEAU_COLORS["tab:blue"] + "55"
    ),
    "query terms": dict(zorder=2, lw=1, color="tab:brown", fc="none"),
    "approximate matches index terms": dict(
        zorder=1,
        color="tab:gray",
        fc="none",
        lw=0.75,
        alpha=0.33,
    ),
    "observations index terms": dict(zorder=1, color="none", fc="tab:pink", alpha=0.33),
    "all observations": dict(
        zorder=0, lw=0.75, color="tab:purple", fc="none", alpha=0.33
    ),
    "matched observations": dict(
        zorder=0, lw=0.75, color="tab:red", fc="none", alpha=0.33
    ),
}

with open(sys.argv[1]) as inf:
    data = json.load(inf)


def wrap(coords):
    i = coords.T[0] < 0
    coords.T[0, i] += 360
    return coords


def observations_polygons(ax):
    polygons = wrap(np.array(list(data["observations"]["polygons"].values())))
    collection = PolyCollection(polygons, closed=True, **styles["all observations"])
    ax.add_collection(collection)
    return polygons


def observations_terms(ax):
    polygons = wrap(np.array(list(data["observations"]["terms"].values())))
    collection = PolyCollection(
        polygons, closed=True, **styles["approximate matches index terms"]
    )
    ax.add_collection(collection)
    return polygons


def matches(ax):
    obsids = data["matches"]
    polygons = wrap(
        np.array([data["observations"]["polygons"][str(obsid)] for obsid in obsids])
    )
    collection = PolyCollection(polygons, closed=True, **styles["matched observations"])
    ax.add_collection(collection)
    return polygons


def ephemeris_polygons(ax):
    polygons = wrap(np.array(list(data["ephemeris"]["polygons"])))
    collection = PolyCollection(polygons, closed=True, **styles["query area"])
    ax.add_collection(collection)
    return polygons


def ephemeris_segments(ax):
    for segment in data["ephemeris"]["segments"]:
        coords = wrap(np.array([(p["ra"], p["dec"]) for p in segment]))
        ax.plot(*coords.T, lw=2, zorder=100)


def ephemeris_terms(ax):
    polygons = wrap(np.array(list(data["ephemeris"]["terms"].values())))
    collection = PolyCollection(polygons, closed=True, **styles["query terms"])
    ax.add_collection(collection)
    return polygons


# def ellipses():
#     phi = np.linspace(0, 360, 361) * u.deg
#     for segment in data["ephemeris segments"]:
#         for ra, dec, a, b, theta in segment:
#             rho = a * b / np.hypot(b * np.cos(phi), a * np.sin(phi)) * u.arcsec
#             yield offset_by(ra * u.deg, dec * u.deg, phi + theta * u.deg, rho)


fig, ax = plt.subplots(num=1, figsize=(10, 10), clear=True)

polygons = observations_polygons(ax)
x = minmax(polygons[..., 0].ravel())
y = minmax(polygons[..., 1].ravel())

# polygons = observations_terms(ax)
# x = np.r_[x, minmax(polygons[..., 0].ravel())]
# y = np.r_[y, minmax(polygons[..., 1].ravel())]

matches(ax)

ephemeris_polygons(ax)
ephemeris_segments(ax)

# labels = "ABCD"
# for k, polygon in enumerate(polygons()):
#     ax.add_patch(polygon)
#     ax.scatter(*polygon.xy.T)

#     for i, j in ([0, 1], [1, 2], [2, 3], [3, 0]):
#         xx, yy = polygon.xy[[i, j]].T
#         ax.plot(xx, yy)
#         ax.text(xx.mean(), yy.mean(), labels[i])

# for ellipse in ellipses():
#     ax.plot(np.r_[ellipse[0], ellipse[0][-1]], np.r_[ellipse[1], ellipse[1][-1]])
#     x.extend(minmax(ellipse[0].value))
#     y.extend(minmax(ellipse[1].value))

# for ra, dec in ephemeris():
#     x.extend(minmax(ra))
#     y.extend(minmax(dec))
#     ax.plot(ra, dec, color="tab:gray", zorder=100)

# ax.add_collection(cells("query terms"))
# ax.add_collection(cells("approximate matches index terms"))
# ax.add_collection(cells("matches index terms"))

# ax.scatter(*test_points.T, marker="o", color="tab:red")

ax.set_xlim((np.max(x) + 1, np.min(x) - 1))
ax.set_ylim((np.min(y) - 1, np.max(y) + 1))
ax.set_aspect(1, "datalim")
ax.minorticks_on()
plt.setp(ax, xlabel="RA (deg)", ylabel="Dec (deg)")
plt.tight_layout()
plt.savefig("query-info.pdf")
