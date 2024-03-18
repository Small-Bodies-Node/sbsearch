# Licensed with the 3-clause BSD license.  See LICENSE for details.

import argparse
from time import monotonic
from collections import defaultdict
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TABLEAU_COLORS
from matplotlib.collections import PatchCollection, PathCollection
from matplotlib.patches import Rectangle
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import SphericalCircle
import astropy.units as u

from . import plot_terms, plot_observations
from ..sbsearch import SBSearch, IntersectionType
from .. import model
from ..target import FixedTarget


# plot styles
styles: Dict[str, dict] = {
    "fixed point": dict(
        zorder=4, marker="*", color="y", edgecolor="k", linewidths=0.5, s=50
    ),
    "fixed area": dict(
        zorder=3, lw=1, color="tab:blue", fc=TABLEAU_COLORS["tab:blue"] + "55"
    ),
    "query cells": dict(zorder=2, lw=1, color="tab:brown", fc="none"),
    "matched cells": dict(zorder=1, color="none", fc="tab:pink", alpha=0.33),
    "observations": dict(zorder=0, lw=0.75, color="tab:red", fc="none", alpha=0.33),
}


def plot_fixed(
    sbs: SBSearch,
    target: FixedTarget,
    sources: List[model.Observation] = [model.Observation],
    ax: Optional[mpl.axes.Axes] = None,
) -> None:
    """Visualize a fixed-target search and results.


    Parameters
    ----------
    sbs : `SBSearch`
        The `SBSearch` object to use for the search.

    target : `FixedTarget`
        The target to search for.

    sources : list of `model.Observation`, optional
        The sources to find in the database.  Set to ``[Observation]`` to search
        all sources.

    ax: `astropy.visualization.wcsaxes.WCSAxes`, optional
        Plot to this axis.

    """

    if ax is None:
        # WCS for plot projection
        wcs = WCS(naxis=2)
        wcs.wcs.crpix = [0, 0]
        wcs.wcs.cdelt = np.array([-1, 1]) / 3600
        wcs.wcs.crval = [target.ra.deg, target.dec.deg]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        wcs.wcs.radesys = "ICRS"
        ax = plt.axes(projection=wcs)

    timestamps: Tuple[str, float] = []

    source: model.Observation
    for sbs.source in sources:
        timestamps.append((f"Searching {sbs.source.__tablename__}", monotonic()))

        # get query terms
        query_terms: Set[str]
        if sbs.intersection_type == IntersectionType.ImageContainsCenter:
            query_terms = set(sbs.indexer.query_point(target.ra.rad, target.dec.rad))
        else:
            query_terms = set(
                sbs.indexer.query_cap(
                    target.ra.rad, target.dec.rad, np.radians(sbs.padding / 60)
                )
            )
        timestamps.append(("* Got query terms", monotonic()))

        # get matching observations from database
        observations: List[model.Observation] = sbs.find_observations_intersecting_cap(
            target
        )
        obs_terms: Set[str] = set(sum([obs.spatial_terms for obs in observations], []))

        timestamps.append(("* Got observations", monotonic()))

        legend_handles: list = []

        if sbs.intersection_type != IntersectionType.ImageContainsCenter:
            ax.add_patch(
                SphericalCircle(
                    target.coordinates,
                    radius=sbs.padding * u.arcmin,
                    transform=ax.get_transform("icrs"),
                    **styles["fixed area"],
                )
            )

        legend_handles.append(
            ax.scatter(
                target.ra.deg,
                target.dec.deg,
                transform=ax.get_transform("icrs"),
                **styles["fixed point"],
            )
        )

        legend_handles.extend(
            [
                plot_terms(
                    ax,
                    query_terms,
                    label="Query",
                    **styles["query cells"],
                ),
                plot_terms(
                    ax,
                    query_terms.intersection(obs_terms),
                    label="Matched cells",
                    **styles["matched cells"],
                ),
                plot_observations(
                    ax,
                    observations,
                    label="Observations",
                    **styles["observations"],
                ),
            ]
        )

        cdelt = wcs.proj_plane_pixel_scales()
        width = max(sbs.padding / 60 * 6, 1) * u.deg
        xlim = (np.r_[-1, 1] * width / 2 / cdelt[0]).to_value("") + wcs.wcs.crpix[0]
        ylim = (np.r_[-1, 1] * width / 2 / cdelt[1]).to_value("") + wcs.wcs.crpix[1]

        plt.setp(
            ax,
            xlim=xlim,
            ylim=ylim,
            xlabel="RA",
            ylabel="Dec",
            aspect="equal",
            adjustable="datalim",
        )

        for a in ax.coords:
            a.display_minor_ticks(True)

        ax.legend(handles=legend_handles, loc="upper left", fontsize="medium")
        plt.tight_layout(pad=0.5)

        timestamps.append(("* Generated plot", monotonic()))

        print()
        t0 = timestamps[0][1]
        last = t0
        for timestamp in timestamps:
            print(
                f"{timestamp[0]:30s} {timestamp[1] - t0:.3f} {timestamp[1] - last:.3f}"
            )
            last = timestamp[1]
