"""
plot S2 text formatted polygons
"""

import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PolyCollection
import astropy.units as u
from astropy.coordinates import offset_by
from mskpy import minmax

with open(sys.argv[1]) as inf:
    data = json.load(inf)

# test_points = []
# with open("test.dat") as inf:
#     for line in inf.read(-1).splitlines():
#         if not line.startswith("["):
#             continue

#         test_points.append(np.array(line[1:-1].split(",")[::-1], float))
# test_points = np.array(test_points)


def polygons():
    # for polygon in np.array(data["query polygons"]):
    #     yield Polygon(polygon, closed=False, color="tab:blue", alpha=0.2)
    collection = PolyCollection(
        np.array(data["query polygons"]), closed=False, color="tab:blue", alpha=0.2
    )
    return collection


def ellipses():
    phi = np.linspace(0, 360, 361) * u.deg
    for segment in data["ephemeris segments"]:
        for ra, dec, a, b, theta in segment:
            rho = a * b / np.hypot(b * np.cos(phi), a * np.sin(phi)) * u.arcsec
            yield offset_by(ra * u.deg, dec * u.deg, phi + theta * u.deg, rho)


def ephemeris():
    for segment in data["ephemeris segments"]:
        yield np.array(segment)[:, :2].T


fig, ax = plt.subplots(num=1, figsize=(10, 10), clear=True)

ax.add_collection(polygons())

# labels = "ABCD"
# for k, polygon in enumerate(polygons()):
#     ax.add_patch(polygon)
#     ax.scatter(*polygon.xy.T)

#     for i, j in ([0, 1], [1, 2], [2, 3], [3, 0]):
#         xx, yy = polygon.xy[[i, j]].T
#         ax.plot(xx, yy)
#         ax.text(xx.mean(), yy.mean(), labels[i])

x = []
y = []
# for ellipse in ellipses():
#     ax.plot(np.r_[ellipse[0], ellipse[0][-1]], np.r_[ellipse[1], ellipse[1][-1]])
#     x.extend(minmax(ellipse[0].value))
#     y.extend(minmax(ellipse[1].value))

for ra, dec in ephemeris():
    x.extend(minmax(ra))
    y.extend(minmax(dec))
    ax.plot(ra, dec, color="tab:gray", zorder=100)

# ax.scatter(*test_points.T, marker="o", color="tab:red")

ax.set_xlim((np.max(x), np.min(x)))
ax.set_ylim((np.min(y), np.max(y)))
ax.set_aspect(1, "datalim")
ax.minorticks_on()
plt.setp(ax, xlabel="RA (deg)", ylabel="Dec (deg)")
plt.tight_layout()
plt.savefig("query-info.pdf")
