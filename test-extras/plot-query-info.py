"""
plot S2 text formatted polygons
"""

import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import TABLEAU_COLORS
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord
import ligo.skymap.plot
import s2geometry as s2
from mskpy import minmax

parser = argparse.ArgumentParser()
parser.add_argument("info_file")
parser.add_argument("--projection", default="mollweide")
# parser.add_argument("--all-sky", action="store_true")
parser.add_argument("--center", type=SkyCoord)
args = parser.parse_args()

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

with open(args.info_file) as inf:
    data = json.load(inf)


def wrap(coords):
    i = coords.T[0] < 0
    coords.T[0, i] += 360
    return coords


def split(polygons, ra_center):
    new_polygons = []
    for poly in polygons:
        ra, dec = poly.T
        if any(ra * u.deg > ra_center + 180 * u.deg) and any(
            ra < ra_center + 180 * u.deg
        ):
            breakpoint()


def observations_polygons(**kwargs):
    polygons = split(
        wrap(np.array(list(data["observations"]["polygons"].values()))), args.center.ra
    )
    collection = PolyCollection(
        polygons, closed=True, **styles["all observations"], **kwargs
    )
    return polygons, collection


def observations_terms(**kwargs):
    polygons = wrap(np.array(list(data["observations"]["terms"].values())))
    collection = PolyCollection(
        polygons, closed=True, **styles["approximate matches index terms"], **kwargs
    )
    return polygons, collection


def matches(**kwargs):
    obsids = data["matches"]
    if len(obsids) == 0:
        return []
    polygons = wrap(
        np.array([data["observations"]["polygons"][str(obsid)] for obsid in obsids])
    )
    collection = PolyCollection(
        polygons, closed=True, **styles["matched observations"], **kwargs
    )
    return polygons, collection


def ephemeris_polygons(**kwargs):
    polygons = wrap(np.array(list(data["ephemeris"]["polygons"])))
    collection = PolyCollection(polygons, closed=True, **styles["query area"], **kwargs)
    return polygons, collection


def ephemeris_segments():
    coords = []
    for segment in data["ephemeris"]["segments"]:
        coords.append(wrap(np.array([(p["ra"], p["dec"]) for p in segment])).T)
    return coords


def ephemeris_terms(**kwargs):
    polygons = wrap(np.array(list(data["ephemeris"]["terms"].values())))
    collection = PolyCollection(
        polygons, closed=True, **styles["query terms"], **kwargs
    )
    return polygons, collection


# def ellipses():
#     phi = np.linspace(0, 360, 361) * u.deg
#     for segment in data["ephemeris segments"]:
#         for ra, dec, a, b, theta in segment:
#             rho = a * b / np.hypot(b * np.cos(phi), a * np.sin(phi)) * u.arcsec
#             yield offset_by(ra * u.deg, dec * u.deg, phi + theta * u.deg, rho)

polygons = wrap(np.array(list(data["observations"]["polygons"].values())))
x = minmax(polygons[..., 0].ravel())
y = minmax(polygons[..., 1].ravel())

fig = plt.figure(1, figsize=(10, 6))
fig.clear()

# WCS for plot projections
wcs = WCS(naxis=2)
# wcs.wcs.
# wcs.wcs.crpix = [0, 0]
# wcs.wcs.cdelt = np.array([-1, 1])
# wcs.wcs.crval = [np.mean((np.min(x), np.max(x))), np.mean((np.min(y), np.max(y)))]
# wcs.wcs.ctype = [f"RA---{args.projection}", f"DEC--{args.projection}"]
# wcs.wcs.radesys = "ICRS"
wcs = WCS(
    {
        "naxis": 2,
        "naxis1": 360,
        "naxis2": 180,
        "crpix1": 180,
        "crpix2": 90,
        "crval1": 180,
        "crval2": 0,
        "cdelt1": -1,
        "cdelt2": 1,
        "ctype1": "RA---AIT",
        "ctype2": "DEC--AIT",
    }
)

# ax = plt.axes(projection=wcs, frame_class=EllipticalFrame)
ax = plt.axes(projection="astro hours {}".format(args.projection), center=args.center)
transform = {"transform": ax.get_transform("world")}

collections = [
    observations_polygons(**transform)[1],
    observations_terms(**transform)[1],
    matches(**transform)[1],
    ephemeris_polygons(**transform)[1],
    ephemeris_terms(**transform)[1],
]
for collection in collections:
    ax.add_collection(collection)

for coords in ephemeris_segments():
    ax.plot(*coords, lw=2, zorder=100, **transform)

# for ellipse in ellipses():
#     ax.plot(np.r_[ellipse[0], ellipse[0][-1]], np.r_[ellipse[1], ellipse[1][-1]])
#     x.extend(minmax(ellipse[0].value))
#     y.extend(minmax(ellipse[1].value))

# if args.all_sky:
#     # ax.set_xlim(-0.5, 360 + 0.5)
#     # ax.set_ylim(-0.5, 180 + 0.5)
#     # ax.set_aspect(1.0)
#     pass
# else:
#     ax.set_xlim((np.max(x) + 1, np.min(x) - 1))
#     ax.set_ylim((np.min(y) - 1, np.max(y) + 1))
#     ax.set_aspect(1, "datalim")

ax.grid()
plt.setp(ax, xlabel="RA (hour angle)", ylabel="Dec (deg)")
plt.tight_layout()
plt.savefig("query-info.pdf")
