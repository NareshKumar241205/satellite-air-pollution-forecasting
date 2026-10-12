"""Report figures from the saved result tables (no model is refitted here)."""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

BLUE, ORANGE, GREY = "#4C78A8", "#F58518", "0.45"


def skill_vs_lead(res, out) -> None:
    v = pd.read_csv(res / "C3_forecast" / "val_by_lead.csv")
    t = pd.read_csv(res / "C3_forecast" / "test2025_by_lead.csv")
    fig, ax = plt.subplots(figsize=(7, 3.8))
    ax.axhline(0, color=GREY, lw=0.8)
    ax.plot(v.days, 100 * v.skill_fused, "-o", ms=3, color=BLUE, label="fused, validation 2023-24")
    ax.plot(t.days, 100 * t.skill_fused, "-o", ms=3, color=ORANGE, label="fused, test 2025 (scored once)")
    ax.plot(t.days, 100 * t.skill_damped, "--", color=ORANGE, lw=1, alpha=0.7, label="damped persistence, test 2025")
    ax.set_xlabel("lead (days)")
    ax.set_ylabel("MSE skill vs climatology (%)")
    ax.set_title("NO2 forecast skill by lead, 18 Tamil Nadu sites (0.04° blocks)", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(out / "fig_skill_vs_lead.png", dpi=150)
    plt.close(fig)


def fan(res, out, sites=("Chennai", "Coimbatore", "Madurai")) -> None:
    f = pd.read_csv(res / "C3_forecast" / "fan_2025_from_last_2024_window.csv", parse_dates=["date"])
    fig, axes = plt.subplots(len(sites), 1, figsize=(7, 2.4 * len(sites)), sharex=True)
    for ax, s in zip(axes, sites):
        g = f[f.site == s]
        ax.fill_between(g.date, g.lo10, g.hi90, color=BLUE, alpha=0.2, label="10-90 % band")
        ax.plot(g.date, g.forecast, color=BLUE, label="forecast (origin: last 2024 window)")
        ax.plot(g.date, g.clim, color=GREY, lw=0.8, ls="--", label="climatology")
        obs = g.observed.where(g.obs_coverage >= 0.5)
        ax.plot(g.date, obs, "o", ms=3, color=ORANGE, label="observed (>=50 % of blocks)")
        ax.set_title(f"{s}: site-mean tropospheric NO2", fontsize=9)
        ax.set_ylabel("µmol/m²")
    axes[0].legend(fontsize=7, frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(out / "fig_fan_2025.png", dpi=150)
    plt.close(fig)


def land_panels(res, out) -> None:
    d = pd.DataFrame(json.loads((res / "B_decision.json").read_text())["layers"])
    d = d[d.available.astype(bool)]
    names = {"dw_label": "DW mode label", "dw_prob_trail3": "DW prob., trailing 3-y median",
             "obt": "Open Buildings Temporal"}
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.4))
    y = np.arange(len(d))
    a1.barh(y, d.mde, color=BLUE)
    a1.axvline(0.5, color=ORANGE, ls="--", lw=1, label="gate 0.5 %/pp")
    a1.axvline(1.0, color=GREY, ls=":", lw=1, label="1.0 %/pp")
    a1.set_yticks(y, [names[x] for x in d.layer], fontsize=8)
    a1.set_xlabel("G1 detectable effect (% NO2 per pp)")
    a1.legend(fontsize=7, frameon=False)
    a2.errorbar(d.lead, y, xerr=[d.lead - d.lead_lo, d.lead_hi - d.lead], fmt="o", color=BLUE, label="lead L(t+1)")
    a2.errorbar(d.WL, y + 0.2, xerr=[d.WL - d.WL_lo, d.WL_hi - d.WL], fmt="s", color=ORANGE, ms=4,
                label="neighbour WL")
    a2.axvline(0, color=GREY, lw=0.8)
    a2.set_yticks(y, [""] * len(d))
    a2.set_xlabel("coefficient, 95 % block-bootstrap CI")
    a2.legend(fontsize=7, frameon=False, loc="lower right")
    fig.suptitle("Land-change design checks per built-up layer (2019-2023, raw Jan-Apr NO2)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "fig_land_layers_g1_lead.png", dpi=150)
    plt.close(fig)


def run(cfg: dict) -> None:
    res = cfg["outputs"]
    out = res / "figures"
    out.mkdir(parents=True, exist_ok=True)
    land_panels(res, out)
    if (res / "C3_forecast" / "test2025_by_lead.csv").exists():
        skill_vs_lead(res, out)
        fan(res, out)
