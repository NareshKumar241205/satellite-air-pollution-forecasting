"""Lead-test diagnostics D1-D6 (Land Change/DIAGNOSTICS_PLAN.md). Tuning years 2019-2023 only."""

from __future__ import annotations

import numpy as np
import pandas as pd

from aqf_land.panel import demean, fit

LEADLAG = ["L", "L_lead", "L_lag", "WL"]


def _sy(d: pd.DataFrame) -> pd.DataFrame:
    return d.assign(sy=d.site + ":" + d.year.astype(str))


def lead_coef(p: pd.DataFrame) -> float:
    return float(fit(p, LEADLAG)[1])


def residuals(p: pd.DataFrame) -> pd.DataFrame:
    """Two-way-FE residuals of y and of L_lead from the lead/lag model (FWL partials)."""
    d = p.dropna(subset=["y", *LEADLAG])
    d = _sy(d[d.groupby("pixel").y.transform("size") >= 2])
    Z = demean(d, ["y", *LEADLAG], ("pixel", "sy"))
    y, X = Z[:, 0], Z[:, 1:]
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    others = [0, 2, 3]
    g, *_ = np.linalg.lstsq(X[:, others], X[:, 1], rcond=None)
    return d.assign(e_y=y - X @ b, e_lead=X[:, 1] - X[:, others] @ g)


def d1_covariates(p: pd.DataFrame, cov: pd.DataFrame) -> pd.DataFrame:
    """Within-site-year correlations of each covariate (year t and t+1) with L(t), dL(t->t+1) and the NO2 residual."""
    c = cov.rename(columns={"city": "site"})
    names = [k for k in ("dw_images", "cloud_prob", "aod055", "brightness") if k in c]
    d = p.merge(c[["site", "row", "col", "year", *names]], on=["site", "row", "col", "year"], how="left")
    nxt = c[["site", "row", "col", "year", *names]].assign(year=lambda x: x.year - 1)
    d = d.merge(nxt, on=["site", "row", "col", "year"], how="left", suffixes=("", "_next"))
    d["dL_next"] = d.L_lead - d.L
    d = _sy(d)
    d["y_res"] = np.nan
    ok = d[["y"]].notna().all(1) & (d.groupby("pixel").y.transform("count") >= 2)
    d.loc[ok, "y_res"] = demean(d[ok], ["y"], ("pixel", "sy"))[:, 0]
    rows = []
    for n in names:
        for cv in (n, f"{n}_next"):
            for target in ("L", "dL_next", "y_res", "L_lead"):
                x = d[[cv, target, "sy"]].dropna()
                if len(x) < 50:
                    continue
                a = x[cv] - x.groupby("sy")[cv].transform("mean")
                b = x[target] - x.groupby("sy")[target].transform("mean")
                rows.append({"covariate": cv, "with": target, "r_within_site_year": float(np.corrcoef(a, b)[0, 1]),
                             "n": len(x)})
    return pd.DataFrame(rows)


def d2_level_dependence(p: pd.DataFrame) -> pd.DataFrame:
    """Per year pair: within-site slope of dL on the pixel's NO2 level (100 log) and on its L level."""
    lvl = p.groupby("pixel").y.mean().rename("y_mean")
    rows = []
    for y0 in sorted(p.year.unique())[:-1]:
        a = p[p.year == y0].set_index("pixel")
        b = p[p.year == y0 + 1].set_index("pixel")
        d = pd.DataFrame({"site": a.site, "dL": b.L - a.L, "L0": a.L}).join(lvl).dropna()
        for c in ("dL", "L0", "y_mean"):
            d[c] = d[c] - d.groupby("site")[c].transform("mean")
        k1 = np.polyfit(d.y_mean, d.dL, 1)[0]
        k2 = np.polyfit(d.L0, d.dL, 1)[0]
        rows.append({"pair": f"{y0}->{y0 + 1}", "dL_per_100logNO2": k1 * 100, "dL_on_L0": k2,
                     "dL_sd_within": d.dL.std(), "n": len(d)})
    return pd.DataFrame(rows)


def d3_pixel_trends(p: pd.DataFrame, it: int = 300) -> dict:
    """Lead/lag model with pixel FE + pixel linear trends + site-year FE (alternating projections)."""
    d = p.dropna(subset=["y", *LEADLAG])
    d = _sy(d[d.groupby("pixel").y.transform("size") >= 3]).reset_index(drop=True)
    Z = d[["y", *LEADLAG]].to_numpy(float).copy()
    t = d.year.to_numpy(float)
    pix = pd.factorize(d.pixel)[0]
    sy = pd.factorize(d.sy)[0]
    npx = pix.max() + 1
    cp = np.bincount(pix, minlength=npx)
    tc = t - (np.bincount(pix, t, npx) / cp)[pix]                 # time centred within pixel
    stt = np.bincount(pix, tc * tc, npx)
    for _ in range(it):
        prev = Z.copy()
        for j in range(Z.shape[1]):                               # remove each pixel's own line in time
            z = Z[:, j] - (np.bincount(pix, Z[:, j], npx) / cp)[pix]
            Z[:, j] = z - (np.bincount(pix, tc * z, npx) / stt)[pix] * tc
        cnt = np.bincount(sy)
        Z -= (np.stack([np.bincount(sy, Z[:, j]) for j in range(Z.shape[1])], 1) / cnt[:, None])[sy]
        if np.abs(Z - prev).max() < 1e-8:
            break
    b, *_ = np.linalg.lstsq(Z[:, 1:], Z[:, 0], rcond=None)
    return {"lead_with_pixel_trends": float(b[1]), "L": float(b[0]), "n": len(d),
            "note": "outcomes 2020-2022: 3 periods per pixel, so 1 residual df per pixel after a pixel trend"}


def d4_leave_one_out(p: pd.DataFrame) -> pd.DataFrame:
    full = lead_coef(p)
    rows = [{"dropped": "none", "lead": full}]
    for s in sorted(p.site.unique()):
        rows.append({"dropped": f"site {s}", "lead": lead_coef(p[p.site != s])})
    for y in (2020, 2021, 2022):
        q = p.copy()
        q.loc[q.year == y, "y"] = np.nan
        rows.append({"dropped": f"outcome year {y}", "lead": lead_coef(q)})
    r = pd.DataFrame(rows)
    r["change_vs_full_pct"] = 100 * (r.lead - full) / full
    return r


def d5_drop_edges(p: pd.DataFrame) -> dict:
    mx = p.groupby("site")[["row", "col"]].transform("max")
    inner = p[(p.row > 0) & (p.col > 0) & (p.row < mx.row) & (p.col < mx.col)]
    full, inn = lead_coef(p), lead_coef(inner)
    return {"lead_full": full, "lead_inner": inn, "change_pct": 100 * (inn - full) / full}


def d6_spatial(p: pd.DataFrame) -> dict:
    r = residuals(p)
    key = r.set_index(["site", "year", "row", "col"])
    out = {}
    for name, (a, b) in {"e_y~e_y": ("e_y", "e_y"), "e_lead~e_lead": ("e_lead", "e_lead"),
                         "e_y~e_lead(neighbour)": ("e_y", "e_lead"),
                         "e_y~e_lead(same pixel)": ("e_y", "e_lead")}.items():
        if "same pixel" in name:
            out[name] = float(np.corrcoef(r[a], r[b])[0, 1])
            continue
        xs, ys = [], []
        for dr, dc in ((0, 1), (1, 0)):
            nb = r[["site", "year", "row", "col", b]].assign(row=r.row - dr, col=r.col - dc)
            m = r[["site", "year", "row", "col", a]].merge(nb, on=["site", "year", "row", "col"],
                                                           suffixes=("", "_nb"))
            xs.append(m[a].to_numpy())
            ys.append(m[f"{b}_nb" if a == b else b].to_numpy())
        out[name] = float(np.corrcoef(np.concatenate(xs), np.concatenate(ys))[0, 1])
    del key
    return out
