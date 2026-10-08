"""Shared figure helpers: consistent style, gas maps with industrial hubs, saving."""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from aqf.data import S5PCube  # noqa: E402

log = logging.getLogger(__name__)
GAS_COLORS = {"NO2": "#c0392b", "CO": "#2c3e50", "SO2": "#d68910"}
plt.rcParams.update({"figure.dpi": 110, "savefig.dpi": 150, "axes.grid": True, "grid.alpha": 0.3,
                     "axes.spines.top": False, "axes.spines.right": False, "font.size": 9})


def save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    log.info("figure %s", path)


def scaled(cfg: dict, gas: str, x):
    return x * cfg["display"][gas]["scale"]


def unit(cfg: dict, gas: str) -> str:
    return cfg["display"][gas]["unit"]


def gas_map(ax, cube: S5PCube, field: np.ndarray, cfg: dict, title: str = "", cmap: str = "viridis",
            vmin=None, vmax=None, hubs: bool = True, label: str = ""):
    """Plot one (H, W) field on the S5P grid in lon/lat, with the industrial hubs marked."""
    dx, dy = np.diff(cube.lon[:2])[0], np.diff(cube.lat[:2])[0]
    extent = [cube.lon[0] - dx / 2, cube.lon[-1] + dx / 2, cube.lat[-1] + dy / 2, cube.lat[0] - dy / 2]
    f = np.where(cube.footprint, field, np.nan)
    im = ax.imshow(f, extent=extent, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
    if hubs:
        for name, (la, lo) in cfg["hubs"].items():
            ax.plot(lo, la, "w^", ms=5, mec="k", mew=0.6)
            ax.annotate(name.split("/")[0], (lo, la), xytext=(3, 3), textcoords="offset points",
                        fontsize=6, color="w", path_effects=[_stroke()])
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("lon (°E)")
    ax.set_ylabel("lat (°N)")
    ax.grid(False)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(label, fontsize=7)
    return im


def _stroke():
    from matplotlib import patheffects
    return patheffects.withStroke(linewidth=1.5, foreground="k")
