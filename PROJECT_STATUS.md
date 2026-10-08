# Project status: 8 Oct 2026

**Satellite air-pollution analysis and forecasting over Chennai** (Sentinel-5P NO2/CO/SO2 + Sentinel-2, 2019–2024).
The project was rebuilt from scratch today. Everything below was re-run from a clean state on 8 Oct 2026 and the
numbers were checked against the reports.

| Check | Result |
|---|---|
| `uv run pytest -q` | **14 passed** (~2 s) |
| `uv run aqf-analysis all` | completed, 53 s, 18 figures + tables |
| `uv run aqf-forecast all` | completed, 27 s, F0 + F1 figures and tables |
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

**Phase 2: Forecasting** (F0 + F1 done, [`Forecasting/PLAN.md`](Forecasting/PLAN.md), [`Forecasting/REPORT.md`](Forecasting/REPORT.md))
- Two tasks:
  - **A**: the next 5-day map.
  - **B**: the **150-day outlook** (30 leads).
- F0 forecastability: anomaly autocorrelation, cross-gas and neighbour predictive correlation, and Sentinel-2 land cover vs pollution with a spatial-block permutation test.
- F1 baselines: climatology, persistence and damped anomaly persistence. Metrics are skill, CRPS, interval coverage and block-bootstrap CIs.
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

| Gas | Task A (+5 days): best baseline | Task A area-mean skill | Task B: beats climatology up to | 80 % interval coverage |
|---|---|---|---|---|
| NO2 | damped persistence **+6.6 % [4.3, 9.9]** | +11 % | **30 days** | 0.79 |
| CO | damped persistence **+9.4 % [0.1, 21]** | +16 % | 5 days | 0.79 |
| SO2 | none (climatology) | — | 0 days | 0.79 |

- Simple persistence (repeating the last map) is 39–83 % *worse* than climatology at +5 days.
- **For leads from about 1 to 5 months, the seasonal climatology is the best forecast available.** The 150-day outlook
  will be the climatology plus a short-lived anomaly correction, with calibrated 10/90 % bands.
- Sentinel-2: mean NDVI vs long-term NO2 gives ρ = −0.53 (block-permutation p = 0.026). Built-up share isn't significant
  (p = 0.14; the naive p of 0.0009 is misleading). S2 explains *where* NO2 is high, not *when*.

## 3. Not done yet

| Phase | Content |
|---|---|
| **F2** (next) | LightGBM + Ridge, pooled over pixels. Inputs: anomaly lags, 3×3 neighbours, trailing 30/90-day mean anomalies, season, static S2, validity masks. Quantile versions give the bands |
| F3 | ConvLSTM, only if F2 shows structure the tree models miss. F0 says this is unlikely |
| F4 | Ablations: S2 off/shuffled, own vs multi-gas, no neighbours, no masks, no 2020 Q2 |
| **F5** | Score 2024 **once** with frozen configs. Then refit on 2019–2024 and produce the **Jan–May 2025 outlook**: a 30-band GeoTIFF per gas (mean, 10 %, 90 %), a 30-row CSV, fan charts, and hotspot series for Manali NO2 and Ennore SO2, plus the next-5-day map |

## 4. Risks and limitations
- **Six years is short.** Trend tests have low power, and only one validation year and one test year exist.
- **SO2 is noisy.** 39 % of its pixels are zeros, 21 % of its windows are unusable, and nothing beats climatology. Its outlook will be climatology only.
- **The 2020 lockdown sits in the training years.** It may inflate the slow NO2 persistence signal, so F4 tests excluding it.
- **Validation numbers drive model choice** and are therefore optimistic. Only the single 2024 score is unbiased.
- **The old pipeline exists only in the system trash.** Don't empty the trash if any of it might still be needed.

## 5. Commands (run from `Scratch/`)

```bash
uv sync
uv run pytest -q
uv run aqf-analysis all      # phase 1, about 1 min, results in "Data Analysis/results/"
uv run aqf-forecast all      # phase 2 F0 + F1, about 30 s, results in Forecasting/results/
```
