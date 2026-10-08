"""F0: how forecastable are the anomalies, and which inputs carry information? Train years only.

1. Anomaly autocorrelation, leads 1..30, for pooled pixels and the area mean. This sets how far ahead
   persistence-type information can reach.
2. Cross-gas lead-lag partial correlation: does gas A at t predict gas B at t+L beyond B's own value at t?
3. Neighbour predictive correlation: does the 3x3 neighbourhood mean at t predict the pixel at t+L beyond its own value?
4. Static Sentinel-2 land cover vs the long-term mean of each gas per S5P pixel (Spearman), with a spatial-block
   permutation test. Shuffling whole 3x3 blocks keeps the spatial autocorrelation that makes naive p-values too small.
"""

from __future__ import annotations

import logging
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.ndimage import convolve

from aqf import s2_static
from aqf_forecast.plotting import GAS_COLORS, save
from aqf_forecast.problem import Problem, fit

log = logging.getLogger(__name__)


def _pairs(a: np.ndarray, mask: np.ndarray, L: int):
    """Values at t and t+L where both windows are in `mask` (a: (T, ...))."""
    both = mask[:-L] & mask[L:]
    return a[:-L][both], a[L:][both]


def _corr(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    return (np.corrcoef(x[ok], y[ok])[0, 1], int(ok.sum())) if ok.sum() > 3 else (np.nan, int(ok.sum()))


def _partial(x, y, z):
    """Correlation of x and y after regressing both on z (with intercept), on jointly valid samples."""
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[ok], y[ok], z[ok]
    Z = np.stack([np.ones_like(z), z], 1)
    rx = x - Z @ np.linalg.lstsq(Z, x, rcond=None)[0]
    ry = y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]
    return float(np.corrcoef(rx, ry)[0, 1]), int(ok.sum())


def run(p: Problem) -> dict:
    fc = p.fcfg
    out = fc["outputs"] / "F0_forecastability"
    out.mkdir(parents=True, exist_ok=True)
    train = p.windows("train")
    f = fit(p, train)
    a = np.where(train[:, None, None, None], f.anom, np.nan)       # (T, G, H, W)
    G = len(p.gases)
    max_lag = fc["forecastability"]["max_lag"]

    # 1. autocorrelation
    with warnings.catch_warnings():                                 # all-NaN windows give NaN, as intended
        warnings.simplefilter("ignore", RuntimeWarning)
        area = np.nanmean(a.reshape(p.T, G, -1), 2)                  # (T, G)
    rows = []
    for g, gas in enumerate(p.gases):
        for L in range(1, max_lag + 1):
            x, y = _pairs(a[:, g], train, L)
            r_pix, n_pix = _corr(x.ravel(), y.ravel())
            xa, ya = _pairs(area[:, g], train, L)
            r_area, n_area = _corr(xa, ya)
            rows.append({"gas": gas, "lag": L, "lag_days": 5 * L, "acf_pixel": r_pix, "n_pixel_pairs": n_pix,
                         "acf_area": r_area, "n_area_pairs": n_area, "white_noise_95": 1.96 / np.sqrt(max(n_area, 1))})
    acf = pd.DataFrame(rows)
    acf.to_csv(out / "acf.csv", index=False)
    reach = []
    for gas in p.gases:
        d = acf[acf.gas == gas].reset_index(drop=True)
        sig = (d.acf_area > d.white_noise_95).values
        first_insig = int(np.argmin(sig)) if (~sig).any() else len(sig)
        reach.append({"gas": gas, "acf_pixel_lag1": d.acf_pixel[0], "acf_area_lag1": d.acf_area[0],
                      "area_acf_significant_up_to_lag": first_insig, "days": 5 * first_insig,
                      "pixel_acf_below_0.1_from_lag": int(d.lag[d.acf_pixel < 0.1].min()) if (d.acf_pixel < 0.1).any() else None})
    reach = pd.DataFrame(reach)
    reach.to_csv(out / "acf_reach.csv", index=False)
    log.info("autocorrelation reach\n%s", reach.round(3).to_string(index=False))

    fig, axes = plt.subplots(1, 2, figsize=(12, 3.8), sharey=True)
    for gas in p.gases:
        d = acf[acf.gas == gas]
        axes[0].plot(d.lag_days, d.acf_pixel, "o-", ms=3, color=GAS_COLORS[gas], label=gas)
        axes[1].plot(d.lag_days, d.acf_area, "o-", ms=3, color=GAS_COLORS[gas], label=gas)
    wn = acf.groupby("lag").white_noise_95.mean()
    axes[1].fill_between(acf.lag_days.unique(), -wn.values, wn.values, color="grey", alpha=0.2, label="95 % noise band")
    for ax, t in zip(axes, ("single pixel (pooled)", "area mean")):
        ax.axhline(0, color="k", lw=0.8)
        ax.set_title(f"Anomaly autocorrelation, {t}, train 2019–2022")
        ax.set_xlabel("lead (days)")
    axes[0].set_ylabel("correlation")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    save(fig, out / "fig_acf.png")

    # 2. cross-gas lead-lag partial correlation (pooled pixels)
    rows = []
    for A, ga in enumerate(p.gases):
        for B, gb in enumerate(p.gases):
            for L in range(1, 7):
                both = train[:-L] & train[L:]
                xa = a[:-L, A][both].ravel()
                yb = a[L:, B][both].ravel()
                zb = a[:-L, B][both].ravel()
                r, n = _partial(xa, yb, zb) if A != B else _corr(zb, yb)
                rows.append({"predictor": ga, "target": gb, "lag": L, "partial_r": r, "n": n})
    cross = pd.DataFrame(rows)
    cross.to_csv(out / "cross_gas_partial.csv", index=False)

    # 3. neighbour predictive correlation (3x3 mean excluding centre, valid pixels only)
    k = np.ones((1, 1, 3, 3))
    k[0, 0, 1, 1] = 0
    fin = np.isfinite(a)
    s = convolve(np.where(fin, a, 0), k, mode="constant")
    c = convolve(fin.astype(float), k, mode="constant")
    with np.errstate(invalid="ignore", divide="ignore"):
        nb = np.where(c > 0, s / c, np.nan)
    rows = []
    for g, gas in enumerate(p.gases):
        for L in range(1, 7):
            both = train[:-L] & train[L:]
            r, n = _partial(nb[:-L, g][both].ravel(), a[L:, g][both].ravel(), a[:-L, g][both].ravel())
            rows.append({"gas": gas, "lag": L, "partial_r_neighbours": r, "n": n})
    neigh = pd.DataFrame(rows)
    neigh.to_csv(out / "neighbour_partial.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 3.8))
    piv = cross[cross.lag == 1].pivot(index="predictor", columns="target", values="partial_r").loc[p.gases, p.gases]
    im = axes[0].imshow(piv.values, cmap="RdBu_r", vmin=-0.3, vmax=0.3)
    for (i, j), v in np.ndenumerate(piv.values):
        axes[0].text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=9)
    axes[0].set_xticks(range(G), [f"{x}(t+5d)" for x in p.gases])
    axes[0].set_yticks(range(G), [f"{x}(t)" for x in p.gases])
    axes[0].grid(False)
    axes[0].set_title("Lead-1 correlation: diagonal = own autocorrelation,\noff-diagonal = partial given target's own value")
    plt.colorbar(im, ax=axes[0], fraction=0.046)
    for gas in p.gases:
        d = neigh[neigh.gas == gas]
        axes[1].plot(5 * d.lag, d.partial_r_neighbours, "o-", color=GAS_COLORS[gas], label=gas)
    axes[1].axhline(0, color="k", lw=0.8)
    axes[1].set_title("3×3 neighbours at t → pixel at t+L, partial on own value")
    axes[1].set_xlabel("lead (days)")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    save(fig, out / "fig_cross_gas_neighbours.png")

    # 4. Sentinel-2 static land cover vs long-term mean (train years)
    feats = s2_static.build(p.cfg, p.cube.lat, p.cube.lon)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        lt_mean = np.nanmean(np.where(train[:, None, None, None], p.values, np.nan), 0)   # (G, H, W)
    fp = p.cube.footprint
    land_ok = fp & (feats["land_frac"] >= 0.5)
    rng = np.random.default_rng(fc["seed"])
    B = fc["forecastability"]["permutation_block"]
    rows = []
    for fname in ("builtup_frac", "ndbi_land", "ndvi_land", "b12_land"):
        for g, gas in enumerate(p.gases):
            x_map, y_map = feats[fname], lt_mean[g]
            rho = stats.spearmanr(x_map[land_ok], y_map[land_ok]).statistic
            null = np.array([_block_perm_rho(x_map, y_map, land_ok, fp, B, rng)
                             for _ in range(fc["forecastability"]["permutation_n"])])
            null = null[np.isfinite(null)]
            p_val = (1 + np.sum(np.abs(null) >= abs(rho))) / (1 + len(null))
            rows.append({"s2_feature": fname, "gas": gas, "spearman_rho": rho, "n_pixels": int(land_ok.sum()),
                         "p_block_permutation": p_val, "p_naive": stats.spearmanr(x_map[land_ok], y_map[land_ok]).pvalue})
    s2 = pd.DataFrame(rows)
    s2.to_csv(out / "s2_vs_longterm_mean.csv", index=False)
    log.info("S2 static vs long-term mean\n%s", s2.round(4).to_string(index=False))

    fig, axes = plt.subplots(1, 3, figsize=(14, 3.8))
    for g, gas in enumerate(p.gases):
        x, y = feats["builtup_frac"][land_ok], lt_mean[g][land_ok] * p.cfg["display"][gas]["scale"]
        axes[g].scatter(x, y, s=12, c=GAS_COLORS[gas], alpha=0.7)
        r = s2[(s2.s2_feature == "builtup_frac") & (s2.gas == gas)].iloc[0]
        axes[g].set_title(f"{gas}: ρ = {r.spearman_rho:.2f}, block-permutation p = {r.p_block_permutation:.3f}", fontsize=9)
        axes[g].set_xlabel("built-up share of land (S2 NDBI > 0)")
        axes[g].set_ylabel(f"mean {gas} ({p.cfg['display'][gas]['unit']})")
    fig.suptitle(f"Static Sentinel-2 land cover vs long-term mean per S5P pixel ({int(land_ok.sum())} mostly-land pixels)")
    fig.tight_layout()
    save(fig, out / "fig_s2_vs_mean.png")
    return {"acf": acf, "reach": reach, "cross": cross, "neighbours": neigh, "s2": s2}


def _block_perm_rho(x_map, y_map, valid, fp, B, rng) -> float:
    """Spearman rho after permuting BxB blocks of x over the footprint's bounding box."""
    rows, cols = np.nonzero(fp)
    r0, c0 = rows.min(), cols.min()
    sub_x = x_map[r0:r0 + 12, c0:c0 + 12]
    nb = 12 // B
    blocks = sub_x.reshape(nb, B, nb, B).transpose(0, 2, 1, 3).reshape(nb * nb, B, B)
    perm = blocks[rng.permutation(nb * nb)]
    # random rotation of each block, so a block never lines up with its original position
    perm = np.stack([np.rot90(b, rng.integers(4)) for b in perm])
    xs = perm.reshape(nb, nb, B, B).transpose(0, 2, 1, 3).reshape(12, 12)
    xp = np.full_like(x_map, np.nan)
    xp[r0:r0 + 12, c0:c0 + 12] = xs
    ok = valid & np.isfinite(xp)
    if ok.sum() < 10:
        return np.nan
    return stats.spearmanr(xp[ok], y_map[ok]).statistic
