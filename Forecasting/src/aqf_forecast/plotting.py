"""Figure helpers for the forecasting phase."""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

log = logging.getLogger(__name__)
GAS_COLORS = {"NO2": "#c0392b", "CO": "#2c3e50", "SO2": "#d68910"}
MODEL_STYLE = {
    "climatology": {"color": "#7f8c8d", "ls": "--", "lw": 1.4, "label": "climatology"},
    "persistence": {"color": "#e67e22", "ls": "-", "lw": 1.4, "label": "persistence"},
    "damped": {"color": "#2471a3", "ls": "-", "lw": 1.8, "label": "damped anomaly persistence"},
}
plt.rcParams.update({"figure.dpi": 110, "savefig.dpi": 150, "axes.grid": True, "grid.alpha": 0.3,
                     "axes.spines.top": False, "axes.spines.right": False, "font.size": 9})


def save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    log.info("figure %s", path)


def scaled(cfg: dict, gas: str, x):
    return x * cfg["display"][gas]["scale"]


def unit(cfg: dict, gas: str) -> str:
    return cfg["display"][gas]["unit"]
