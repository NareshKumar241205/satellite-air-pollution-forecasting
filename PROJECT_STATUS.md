# Project status: 8 Oct 2026

**Satellite air-pollution analysis and forecasting over Chennai** (Sentinel-5P NO2/CO/SO2 + Sentinel-2, 2019–2024).
The project was rebuilt from scratch today. Everything below was re-run from a clean state on 8 Oct 2026 and the
numbers were checked against the reports.

| Check | Result |
|---|---|
| `uv run pytest -q` | **16 passed** (~3 s) |
| `uv run aqf-analysis all` | completed, 53 s, 18 figures + tables |
| `uv run aqf-forecast all` | completed: F0 + F1 about 30 s, F2 learned models about 3 min |
| Git | `main`, local, working tree clean before this file |

## 1. What has been done

**Clean-up and rebuild**
- Inventoried the old `Scratch/`. At the user's request it was reset to only `Dataset/`. The old code is in the system trash and can be restored.
- Removed duplicate files (byte-identical `test_*_20200301.tif`), orphaned run folders and caches. Nothing was hard-deleted; all of it went through `trash-put`.
- Rebuilt the project as a uv workspace with one folder per phase:
  - `src/aqf/`: shared core (data loading, quality rules, climatology, S2 features)
  - `Data Analysis/`: phase 1
  - `Forecasting/`: phase 2
- Git history: 6 commits, one per step.

**Phase 1: Data Analysis** (done, [`Data Analysis/REPORT.md`](<Data Analysis/REPORT.md>))
- Inventory and integrity: 438 five-day windows, one grid per sensor, SHA-256 provenance.
- Quality: **exact zeros are missing retrievals, not clean air**, so they're treated as missing.
- Trends: STL decomposition, seasonal Mann–Kendall / Sen slope, per-pixel trend maps with FDR correction.
- COVID lockdown effect, spatial hotspots, industrial hubs, and a Sentinel-2 summary.

**Phase 2: Forecasting** (F0, F1, F2 done, [`Forecasting/PLAN.md`](Forecasting/PLAN.md), [`Forecasting/REPORT.md`](Forecasting/REPORT.md))
- Two tasks:
  - **A**: the next 5-day map.
  - **B**: the **150-day outlook** (30 leads).
- F0 forecastability: anomaly autocorrelation, cross-gas and neighbour predictive correlation, and Sentinel-2 land cover vs pollution with a spatial-block permutation test.
- F1 baselines: climatology, persistence and damped anomaly persistence. Metrics are skill, CRPS, interval coverage and block-bootstrap CIs.
- F2 learned models: LightGBM and Ridge per gas, pooled over pixels. Task A is a lead-1 model and task B a multi-horizon model with the lead as a feature. 27 features, including 30/60/90-day trailing anomalies. A no-leakage test covers the features, and there's a 2020 Q2 (lockdown) robustness check.
- Fitted on train 2019–2022 and scored on validation 2023. **The 2024 test year is untouched.**

## 2. Key results

### Data (phase 1)
| | NO2 | CO | SO2 |
|---|---|---|---|
| Usable 5-day windows (of 438) | 411 (94 %) | 409 (93 %) | 344 (79 %) |
| Seasonal share of variance (STL) | 0.55 | **0.83** | 0.44 |
| Trend 2019–2024 (seasonal MK) | +2.0 %/yr, p = 0.12 | +0.4 %/yr, p = 0.36 | −0.5 %/yr, p = 0.79 |
| 2020 lockdown vs other years | **−35.5 %** | −0.7 % | −21 % |

- No gas has a statistically significant trend over 2019–2024.
- The NO2 hotspot is on the Manali–Ambattur axis. The SO2 maximum is in the far north-east, toward Ennore. CO is spatially flat (±2 %).

### Forecasting, validation 2023 (phase 2)
Skill means MSE skill against the per-pixel seasonal climatology. Positive means better than the seasonal average. The brackets are 95 % CIs.

**Task A: next 5-day map**

| Gas | Best baseline (damped persistence) | Best learned model (F2) | Learned vs damped | Area-mean skill (learned) |
|---|---|---|---|---|
| NO2 | +6.6 % [4.3, 9.9] | **Ridge +12.7 % [3.6, 22.9]** | +6.5 % [−1.3, 15.0] (n.s.) | **+35 %** |
| CO | +9.4 % [0.1, 21] | **LightGBM +15.5 % [2.9, 29.3]** | +6.8 % [−1.6, 14.4] (n.s.) | **+23 %** |
| SO2 | −0.3 % | +0.2 % | negligible | +1 % |

**Task B: 150-day outlook**

| Gas | Beats climatology up to (baselines) | Beats climatology up to (learned) | Mean skill over 5–150 d |
|---|---|---|---|
| NO2 | 30 days | **130 days (Ridge, also beats damped at leads 2–26); 150 days (LightGBM)** | Ridge +5.9 % |
| CO | 5 days | 10 days (Ridge) | +1.2 % |
| SO2 | 0 | 0 | ≈ 0 |

- The NO2 long-lead skill comes from slow, area-wide anomalies (60/90-day trailing means). It **survives
  excluding the 2020 lockdown** from training (Ridge +5.6 % → +6.9 % at 35–150 days).
- The 80 % intervals from out-of-sample residuals are close to calibrated (coverage 0.76–0.83). Quantile-LightGBM bands under-cover (0.73–0.77).
- ConvLSTM (F3) isn't justified: the signal is area-wide, slow and near-linear.
- Sentinel-2: mean NDVI vs long-term NO2 gives ρ = −0.53 (block-permutation p = 0.026). The S2 features aren't among the top model features.

## 3. Not done yet

| Phase | Content |
|---|---|
| **F4** (next) | Ablations: S2 off/shuffled, own vs multi-gas, no neighbours, no masks (the 2020 Q2 check is done in F2) |
| **F5** | Score 2024 **once** with frozen configs. Candidate models: NO2 Ridge (A and B), CO LightGBM (A) / climatology (B), SO2 climatology. Then refit on 2019–2024 and produce the **Jan–May 2025 outlook**: a 30-band GeoTIFF per gas (mean, 10 %, 90 %), a 30-row CSV, fan charts, and hotspot series for Manali NO2 and Ennore SO2, plus the next-5-day map |

## 4. Risks and limitations
- **Six years is short.** Trend tests have low power, and only one validation year and one test year exist.
- **SO2 is noisy.** 39 % of its pixels are zeros, 21 % of its windows are unusable, and nothing beats climatology. Its outlook will be climatology only.
- **2023 (validation) was the highest-NO2 year**, which favours 'the anomaly persists' models. 2024 fell 15 %, so the single test will show whether the long-lead NO2 skill holds.
- **Learned models are worse than climatology for NO2 in the NE monsoon (Oct–Dec).**
- **Validation numbers drive model choice** and are therefore optimistic. Only the single 2024 score is unbiased.
- **The old pipeline exists only in the system trash.** Don't empty the trash if any of it might still be needed.

## 5. Commands (run from `Scratch/`)

```bash
uv sync
uv run pytest -q
uv run aqf-analysis all      # phase 1, about 1 min, results in "Data Analysis/results/"
uv run aqf-forecast all      # phase 2 F0 + F1 + F2, about 3.5 min, results in Forecasting/results/
```
