"""Gate G1, the power check: is there enough land change across the 12 x 12 pixels to detect an effect on NO2?

Only tuning years are used (configs/land.yaml). No land-NO2 slope is estimated here: G1 asks whether the design
can see an effect before looking for one.

For each city and for all cities pooled:
1. The distribution of dL (built-up fraction change, pp) over pixel-year pairs, and the share of pairs whose
   |dL| exceeds `dL_meaningful_pp`.
2. The noise the effect must stand out from: the residual of dNO2 after removing each city-year mean (the
   fixed effects absorb anything common to a city and year: weather, lockdown, retrieval-version steps).
3. The spread of dL after the same fixed effects, which is the only variation that identifies the slope.
4. The effective sample size, discounting spatial autocorrelation of the dNO2 residual between neighbouring
   pixels: N_eff = N (1 - rho) / (1 + rho), with rho the lag-1 neighbour correlation.
5. dL measurement noise. Built-up loss is rare in reality, so the negative side of dL is treated as noise:
   noise variance = mean(dL^2 | dL < 0), assuming symmetric noise and true dL >= 0. The signal spread is
   sd_signal = sqrt(max(var(dL within) - noise variance, 0)). Classification noise inflates sd(dL) while
   attenuating the slope, so the raw spread must not be used for power.
6. A phenology check: the within-city-year correlation of dL with the vegetation change. Near -1 means dL is
   greenness flipping with rainfall, not land use.
7. The minimum detectable slope MDE = (z_{1-alpha/2} + z_power) * sd(e) / (sd_signal * sqrt(N_eff)), in % NO2
   per pp, and the true-dL spread the design would need for the assumed effect:
   sd_required = (z_{1-alpha/2} + z_power) * sd(e) / (effect * sqrt(N_eff)).
G1 passes when the pooled MDE is at most the assumed effect `effect_pct_per_pp`.
Sensitivity only (row 'ALL_weighted', never the gate): cities combined by precision, MDE = z / sqrt(sum_c N_eff_c
sd_signal_c^2 / sd(e)_c^2), which lets growing cities count more than static ones with their own noise level.
Everything is computed for consecutive-year pairs and for one long difference over the tuning years.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import norm


def neighbour_rho(df: pd.DataFrame, col: str) -> float:
    """Correlation of `col` between horizontally and vertically adjacent pixels of the same city and year pair."""
    a, b = [], []
    for _, g in df.groupby(["city", "y0"]):
        grid = g.pivot(index="row", columns="col", values=col).to_numpy()
        for x, y in ((grid[:, :-1], grid[:, 1:]), (grid[:-1], grid[1:])):
            ok = np.isfinite(x) & np.isfinite(y)
            a.append(x[ok])
            b.append(y[ok])
    a, b = np.concatenate(a), np.concatenate(b)
    return float(np.corrcoef(a, b)[0, 1]) if len(a) > 2 else np.nan


def power_row(df: pd.DataFrame, cfg: dict) -> dict:
    g1 = cfg["g1"]
    d = df.dropna(subset=["dL", "dNO2"]).copy()
    if len(d) < 10:
        return {"n": len(d)}
    fe = d.groupby(["city", "y0"])
    d["e"] = d.dNO2 - fe.dNO2.transform("mean")
    d["dL_w"] = d.dL - fe.dL.transform("mean")
    rho = neighbour_rho(d, "e")
    n_eff = len(d) * (1 - rho) / (1 + rho) if np.isfinite(rho) and rho > 0 else len(d)
    z = norm.ppf(1 - g1["alpha"] / 2) + norm.ppf(g1["power"])
    sd_e, sd_l = d.e.std(ddof=1), d.dL_w.std(ddof=1)
    neg = d.dL[d.dL < 0]
    noise_var = float((neg ** 2).mean()) if len(neg) else 0.0
    sd_signal = float(np.sqrt(max(sd_l ** 2 - noise_var, 0.0)))
    mde = z * sd_e / (sd_signal * np.sqrt(n_eff)) if sd_signal > 0 else np.inf
    q = d.dL.quantile([0.05, 0.25, 0.5, 0.75, 0.95]).to_numpy()
    veg = d.d_vegetation - fe.d_vegetation.transform("mean") if "d_vegetation" in d else None
    return {
        "n": len(d), "n_pixels": d.groupby(["city", "row", "col"]).ngroups,
        "dL_p05": q[0], "dL_p25": q[1], "dL_median": q[2], "dL_p75": q[3], "dL_p95": q[4],
        "dL_sd": d.dL.std(ddof=1), "dL_sd_within": sd_l,
        "share_changed": float((d.dL.abs() > g1["dL_meaningful_pp"]).mean()),
        "share_gain": float((d.dL > g1["dL_meaningful_pp"]).mean()),
        "share_loss": float((d.dL < -g1["dL_meaningful_pp"]).mean()),
        "dL_noise_sd": np.sqrt(noise_var), "dL_signal_sd": sd_signal,
        "corr_dL_dveg": float(np.corrcoef(d.dL_w, veg)[0, 1]) if veg is not None else np.nan,
        "dNO2_resid_sd": sd_e, "rho_neighbour": rho, "n_eff": n_eff,
        "mde_pct_per_pp": mde, "effect_pct_per_pp": g1["effect_pct_per_pp"],
        "sd_required_pp": z * sd_e / (g1["effect_pct_per_pp"] * np.sqrt(n_eff)),
        "passes": bool(mde <= g1["effect_pct_per_pp"]),
    }


def long_difference(p: pd.DataFrame) -> pd.DataFrame:
    """One row per pixel: first -> last tuning year, the sum of the consecutive-year changes (exact for logs)."""
    n = p.y0.nunique()
    cols = ["dL", "dNO2"] + [c for c in p if c.startswith("d_")]
    g = p.groupby(["city", "row", "col"])[cols].sum(min_count=n).reset_index()
    g["y0"], g["y1"] = p.y0.min(), p.y1.max()
    return g


def run_g1(pairs: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Rows: every city plus 'ALL', for four sets: all consecutive tuning pairs, without the lockdown year,
    without the pair that straddles the TROPOMI v1 -> v2 processor change, and one long difference."""
    tune = cfg["years"]["tuning"]
    p = pairs[pairs.y0.isin(tune) & pairs.y1.isin(tune)]
    lock = cfg["years"]["lockdown_year"]
    sets = {
        "all_tuning_pairs": p,
        "no_lockdown_year": p[(p.y0 != lock) & (p.y1 != lock)],
        "no_v1_v2_pair": p[~((p.y0 == 2021) & (p.y1 == 2022))],
        f"long_{tune[0]}_{tune[-1]}": long_difference(p),
    }
    rows = []
    for name, s in sets.items():
        for city, g in s.groupby("city"):
            rows.append({"pairs": name, "city": city, **power_row(g, cfg)})
        rows.append({"pairs": name, "city": "ALL", **power_row(s, cfg)})
        c = pd.DataFrame([r for r in rows if r["pairs"] == name and r["city"] not in ("ALL", "ALL_weighted")])
        info = (c.n_eff * c.dL_signal_sd ** 2 / c.dNO2_resid_sd ** 2).sum()
        z = norm.ppf(1 - cfg["g1"]["alpha"] / 2) + norm.ppf(cfg["g1"]["power"])
        mde = z / np.sqrt(info) if info > 0 else np.inf
        rows.append({"pairs": name, "city": "ALL_weighted", "n": int(c.n.sum()), "mde_pct_per_pp": mde,
                     "effect_pct_per_pp": cfg["g1"]["effect_pct_per_pp"],
                     "passes": bool(mde <= cfg["g1"]["effect_pct_per_pp"])})
    return pd.DataFrame(rows)
