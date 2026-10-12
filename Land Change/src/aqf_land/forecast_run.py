"""Run DECISION_RULES.md C-E: validation ladder, ablations, long range, fusion, then the one-time 2025 test."""

from __future__ import annotations

import json
import logging

import numpy as np
import pandas as pd

from aqf.data import ROOT
from aqf_land import forecast as fc

log = logging.getLogger("aqf_land.forecast")
TRAIN, VAL, TEST = (2019, 2022), (2023, 2024), (2025, 2025)
HMAX = 30


def _in(years: np.ndarray, span: tuple[int, int]) -> np.ndarray:
    return (years >= span[0]) & (years <= span[1])


def load(cfg: dict) -> fc.Data:
    from aqf_land.layers import CLASSES, land_from_csv
    from aqf_land.tnaq import read_cities

    no2 = read_cities(cfg, "NO2")
    era5 = pd.read_csv(cfg["outputs"] / "P1_era5" / "era5_city_windows.csv", parse_dates=["window"])
    k = CLASSES.index("builtup")
    land = {}
    for box in no2["boxes"]:
        shape = (box.rows // cfg["grid"]["block"], box.cols // cfg["grid"]["block"])
        lf = land_from_csv(ROOT / "Land Change/data/ee_landfrac", box.city, shape, list(range(2019, 2026)),
                           cfg["grid"]["min_valid_frac"])
        land[box.city] = {y: v[k] * 100 for y, v in lf.items()}
        land[box.city][2018] = land[box.city][2019]       # documented fallback for origins in Jan-Apr 2019
    return fc.Data(cfg, no2, land, era5)


def h1_rows(a: np.ndarray, years: np.ndarray, span) -> tuple[np.ndarray, np.ndarray]:
    """Origins t whose target t+1 lies in span and is valid: (t, p) index arrays."""
    T = len(a)
    tgt_ok = np.zeros_like(a, bool)
    tgt_ok[:-1] = np.isfinite(a[1:]) & _in(years[1:], span)[:, None]
    return np.nonzero(tgt_ok)


def to_grid(shape, t, p, v) -> np.ndarray:
    out = np.full(shape, np.nan)
    out[t + 1, p] = v
    return out


def fit_ridge(X, y, Xv, yv, alphas=(0.1, 1, 10, 100, 1000)):
    mu, sd = X.mean(0), X.std(0) + 1e-9
    Z, Zv = (X - mu) / sd, (Xv - mu) / sd
    best = None
    for al in alphas:
        A = Z.T @ Z + al * np.eye(Z.shape[1])
        w = np.linalg.solve(A, Z.T @ (y - y.mean()))
        mse = (((Zv @ w + y.mean()) - yv) ** 2).mean()
        if best is None or mse < best[0]:
            best = (mse, al, w, y.mean())
    return best[1], (lambda Q: ((Q - mu) / sd) @ best[2] + best[3])


def fit_lgbm(X, y, Xv, yv, rounds=None):
    import lightgbm as lgb

    dtr = lgb.Dataset(X, y)
    if rounds is None:
        dv = lgb.Dataset(Xv, yv, reference=dtr)
        m = lgb.train(fc.LGB_PARAMS, dtr, 2000, valid_sets=[dv], callbacks=[lgb.early_stopping(50, verbose=False)])
        return m, m.best_iteration
    return lgb.train(fc.LGB_PARAMS, dtr, rounds), rounds


def cnn_try(d, a, years, rng):
    """Small fully-convolutional net on site grids (anomaly lags 0-5 + missing flags + day of year). None if no torch."""
    try:
        import torch
        import torch.nn as nn
    except ImportError:
        return None, "not run: PyTorch is not installed in this environment"
    torch.manual_seed(0)
    grids = []
    for s in range(len(d.sites)):
        m = d.meta[d.site_idx == s]
        H, W = m.row.max() + 1, m.col.max() + 1
        grids.append((np.flatnonzero(d.site_idx == s), H, W))
    T = len(a)

    def tensors(span):
        xs, ys, ms, keys = [], [], [], []
        for s, (idx, H, W) in enumerate(grids):
            g = a[:, idx].reshape(T, H, W)
            for t in range(6, T - 1):
                if not (span[0] <= years[t + 1] <= span[1]):
                    continue
                lags = g[t - 5:t + 1][::-1]
                ch = np.concatenate([np.nan_to_num(lags), np.isnan(lags)], 0)
                tgt = g[t + 1]
                xs.append(np.pad(ch, ((0, 0), (0, 12 - H), (0, 12 - W))))
                ys.append(np.pad(np.nan_to_num(tgt), ((0, 12 - H), (0, 12 - W))))
                ms.append(np.pad(np.isfinite(tgt), ((0, 12 - H), (0, 12 - W))))
                keys.append((s, t))
        return (torch.tensor(np.array(xs), dtype=torch.float32), torch.tensor(np.array(ys), dtype=torch.float32),
                torch.tensor(np.array(ms)), keys)

    Xtr, Ytr, Mtr, _ = tensors(TRAIN)
    Xv, Yv, Mv, kv = tensors(VAL)
    net = nn.Sequential(nn.Conv2d(12, 32, 3, padding=1), nn.ReLU(), nn.Conv2d(32, 32, 3, padding=1), nn.ReLU(),
                        nn.Conv2d(32, 1, 3, padding=1))
    opt = torch.optim.Adam(net.parameters(), 1e-3)
    best, best_state, bad = np.inf, None, 0
    for ep in range(200):
        net.train()
        perm = torch.randperm(len(Xtr))
        for i in range(0, len(perm), 64):
            b = perm[i:i + 64]
            loss = (((net(Xtr[b])[:, 0] - Ytr[b]) ** 2) * Mtr[b]).sum() / Mtr[b].sum()
            opt.zero_grad()
            loss.backward()
            opt.step()
        net.eval()
        with torch.no_grad():
            v = float(((((net(Xv)[:, 0] - Yv) ** 2) * Mv).sum() / Mv.sum()))
        if v < best - 1e-6:
            best, best_state, bad = v, {k: x.clone() for k, x in net.state_dict().items()}, 0
        else:
            bad += 1
            if bad >= 10:
                break
    net.load_state_dict(best_state)
    with torch.no_grad():
        pv = net(Xv)[:, 0].numpy()
    out = np.full(a.shape, np.nan)
    for (s, t), pr in zip(kv, pv):
        idx, H, W = grids[s]
        out[t + 1, idx] = pr[:H, :W].ravel()
    return out, f"run: {ep + 1} epochs, best validation MSE {best:.3f}"


def isotonic_decreasing(g: np.ndarray) -> np.ndarray:
    """Pool-adjacent-violators for a non-increasing sequence."""
    vals = [[v, 1] for v in g]
    i = 0
    while i < len(vals) - 1:
        if vals[i][0] < vals[i + 1][0]:
            v = (vals[i][0] * vals[i][1] + vals[i + 1][0] * vals[i + 1][1]) / (vals[i][1] + vals[i + 1][1])
            vals[i] = [v, vals[i][1] + vals[i + 1][1]]
            del vals[i + 1]
            i = max(i - 1, 0)
        else:
            i += 1
    return np.concatenate([[v] * n for v, n in vals])


def lead_arrays(a, alast, phi, years, span):
    """For each lead h: origins t with target t+h in span; returns dict h -> (t, p, damped anomaly, obs anomaly)."""
    T = len(a)
    out = {}
    for h in range(1, HMAX + 1):
        ok = np.zeros_like(a, bool)
        ok[:T - h] = np.isfinite(a[h:]) & _in(years[h:], span)[:, None]
        t, p = np.nonzero(ok)
        out[h] = (t, p, phi[h] * alast[t, p], a[t + h, p])
    return out


def run(cfg: dict, rescore: bool = False) -> None:
    out = cfg["outputs"] / "C3_forecast"
    out.mkdir(parents=True, exist_ok=True)
    rules = json.loads((cfg["outputs"] / "B_decision.json").read_text())
    if rules["B"] != 0:
        raise NotImplementedError("B != 0: the confirmatory land term is not wired in; see B_decision.json")
    rng = np.random.default_rng(0)
    d = load(cfg)
    years = d.years()
    T, P = d.Y.shape
    tr_w = _in(years, TRAIN)

    # ---------------- validation stage: fit on train only ----------------
    clim = fc.climatology(d, tr_w)
    a = d.Y - clim
    alast = fc.last_valid(a)
    phi = fc.damping(a, alast, tr_w)
    X, names = fc.features(d, a, fc.land_year)
    tt, pt = h1_rows(a, years, TRAIN)
    tv, pv = h1_rows(a, years, VAL)
    Xtr, ytr, Xv, yv = X[tt, pt], a[tt + 1, pt], X[tv, pv], a[tv + 1, pv]
    yv_grid = np.where(_in(years, VAL)[:, None], d.Y, np.nan)
    preds = {"climatology": clim.copy(), "damped persistence": to_grid((T, P), tv, pv, phi[1] * alast[tv, pv]) + clim}
    alpha, ridge = fit_ridge(Xtr, ytr, Xv, yv)
    preds["ridge"] = to_grid((T, P), tv, pv, ridge(Xv)) + clim
    lgbm, best_iter = fit_lgbm(Xtr, ytr, Xv, yv)
    preds["lightgbm"] = to_grid((T, P), tv, pv, lgbm.predict(Xv)) + clim
    cnn, cnn_note = cnn_try(d, a, years, rng)
    if cnn is not None:
        preds["cnn"] = cnn + clim
    mask_v = np.zeros((T, P), bool)
    mask_v[tv + 1, pv] = True
    val_rows = []
    for name, pr in preds.items():
        pr = np.where(mask_v, pr, np.nan)
        val_rows.append({"model": name, **fc.scores(np.where(mask_v, d.Y, np.nan), pr, clim, d.dates, d.site_idx, rng)})
    val = pd.DataFrame(val_rows)
    val.to_csv(out / "val_h1_ladder.csv", index=False)
    chosen = val.loc[val.skill.idxmax(), "model"]
    log.info("validation h=1: %s; chosen %s", val[["model", "skill"]].round(3).to_dict("records"), chosen)

    # ablations (LightGBM, validation)
    abl = []
    for label, kw in {"full": {}, "ERA5 off": {"use_era5": False}, "land off": {"use_land": False},
                      "land shuffled": {"shuffle_land": np.random.default_rng(1)}}.items():
        Xa, _ = (X, names) if not kw else fc.features(d, a, fc.land_year, **kw)
        m, _ = fit_lgbm(Xa[tt, pt], ytr, Xa[tv, pv], yv)
        pr = np.where(mask_v, to_grid((T, P), tv, pv, m.predict(Xa[tv, pv])) + clim, np.nan)
        abl.append({"variant": label, **fc.scores(np.where(mask_v, d.Y, np.nan), pr, clim, d.dates, d.site_idx, rng,
                                                    n_boot=200)})
    pd.DataFrame(abl).to_csv(out / "val_ablations_lightgbm.csv", index=False)

    # long range + fusion weights on validation
    Ahat1 = {"climatology": np.zeros(len(tv)), "damped persistence": phi[1] * alast[tv, pv], "ridge": ridge(Xv),
             "lightgbm": lgbm.predict(Xv)}
    if cnn is not None:
        Ahat1["cnn"] = cnn[tv + 1, pv]
    la = lead_arrays(a, alast, phi, years, VAL)
    g = np.zeros(HMAX + 1)
    lead_rows = []
    for h in range(1, HMAX + 1):
        t, p, A, obs = la[h]
        if h == 1:
            A = Ahat1[chosen]
        den = (A ** 2).sum()
        g[h] = float(np.clip((A * obs).sum() / den, 0, 1)) if den > 0 else 0.0
    g[1:] = isotonic_decreasing(g[1:])
    q = {}
    for h in range(1, HMAX + 1):
        t, p, A, obs = la[h]
        if h == 1:
            A = Ahat1[chosen]
        resid = obs - g[h] * A
        q[h] = (float(np.percentile(resid, 10)), float(np.percentile(resid, 90)))
        dmp = 1 - ((obs - la[h][2]) ** 2).mean() / (obs ** 2).mean()
        fus = 1 - (resid ** 2).mean() / (obs ** 2).mean()
        cover = float(((resid >= q[h][0]) & (resid <= q[h][1])).mean())
        lead_rows.append({"lead": h, "days": 5 * h, "skill_damped": dmp, "skill_fused": fus, "g": g[h],
                          "phi": phi[h], "cover_10_90": cover, "n": len(obs)})
    pd.DataFrame(lead_rows).to_csv(out / "val_by_lead.csv", index=False)
    frozen = {"chosen_h1": chosen, "ridge_alpha": alpha, "lgbm_rounds": int(best_iter), "g": g[1:].tolist(),
              "q10_q90": {h: q[h] for h in q}, "features": names, "cnn": cnn_note, "B": 0,
              "note": "g(h) and intervals are kept from the validation stage for the final refit (in-sample "
                      "shrinkage of a refit ML model would be biased); every other part is refitted on 2019-2024."}
    (out / "FROZEN.json").write_text(json.dumps(frozen, indent=2))

    # ---------------- final: refit on 2019-2024, score 2025 once ----------------
    flag = out / "TEST_SCORED.json"
    if flag.exists() and not rescore:
        log.info("2025 already scored once (%s); not rescoring", flag)
        return
    tv_w = _in(years, (TRAIN[0], VAL[1]))
    clim2 = fc.climatology(d, tv_w)
    a2 = d.Y - clim2
    alast2 = fc.last_valid(a2)
    phi2 = fc.damping(a2, alast2, tv_w)
    X2, _ = fc.features(d, a2, fc.land_year)
    t2, p2 = h1_rows(a2, years, (TRAIN[0], VAL[1]))
    ts, ps = h1_rows(a2, years, TEST)
    if chosen == "lightgbm":
        m2, _ = fit_lgbm(X2[t2, p2], a2[t2 + 1, p2], None, None, rounds=best_iter)
        A1 = m2.predict(X2[ts, ps])
    elif chosen == "ridge":
        _, r2 = fit_ridge(X2[t2, p2], a2[t2 + 1, p2], X2[t2, p2], a2[t2 + 1, p2], alphas=(alpha,))
        A1 = r2(X2[ts, ps])
    elif chosen == "damped persistence":
        A1 = phi2[1] * alast2[ts, ps]
    elif chosen == "climatology":
        A1 = np.zeros(len(ts))
    else:
        raise NotImplementedError("CNN refit not implemented")
    lt = lead_arrays(a2, alast2, phi2, years, TEST)
    rows = []
    for h in range(1, HMAX + 1):
        t, p, A, obs = lt[h]
        if h == 1:
            A = A1
        resid = obs - g[h] * A
        rows.append({"lead": h, "days": 5 * h, "skill_fused": 1 - (resid ** 2).mean() / (obs ** 2).mean(),
                     "skill_damped": 1 - ((obs - lt[h][2]) ** 2).mean() / (obs ** 2).mean(),
                     "rmse_fused": float(np.sqrt((resid ** 2).mean())), "rmse_clim": float(np.sqrt((obs ** 2).mean())),
                     "mae_fused": float(np.abs(resid).mean()),
                     "cover_10_90": float(((resid >= q[h][0]) & (resid <= q[h][1])).mean()), "n": len(obs)})
    test = pd.DataFrame(rows)
    test.to_csv(out / "test2025_by_lead.csv", index=False)
    mask_s = np.zeros((T, P), bool)
    mask_s[ts + 1, ps] = True
    pred1 = to_grid((T, P), ts, ps, g[1] * A1) + clim2
    s1 = fc.scores(np.where(mask_s, d.Y, np.nan), np.where(mask_s, pred1, np.nan), clim2, d.dates, d.site_idx, rng)
    flag.write_text(json.dumps({"scored_once": True, "h1": s1, "chosen_h1": chosen,
                                "mean_skill_leads_2_30": float(test.skill_fused[1:].mean())}, indent=2))
    log.info("2025 test (once): h=1 skill %.3f [%.3f, %.3f]", s1["skill"], s1["skill_lo"], s1["skill_hi"])
