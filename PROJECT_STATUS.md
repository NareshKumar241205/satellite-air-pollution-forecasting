# Project status: 8 Oct 2026

**Satellite air-pollution analysis and forecasting over Chennai** (Sentinel-5P NO2/CO/SO2 + Sentinel-2, 2019–2024).
The project was rebuilt from scratch today. Everything below was re-run from a clean state on 8 Oct 2026 and the
numbers were checked against the reports.

| Check | Result |
|---|---|
| `uv run pytest -q` | **21 passed** (~2 s) |
| `uv run aqf-analysis all` | completed, 53 s, 18 figures + tables |
| `uv run aqf-forecast all` | completed: F0 + F1 about 30 s, F2 about 3 min, F4 about 30 s |
| `uv run aqf-forecast final` | 2024 scored **once** + Jan–May 2025 outlook, about 1 min |
| Git | `main`, local, working tree clean before this file |

## 1. What has been done

**Clean-up and rebuild**
- Inventoried the old `Scratch/`. At the user's request it was reset to only `Dataset/`. The old code is in the system trash and can be restored.
- Removed duplicate files (byte-identical `test_*_20200301.tif`), orphaned run folders and caches. Nothing was hard-deleted; all of it went through `trash-put`.
- Rebuilt the project as a uv workspace with one folder per phase:
  - `src/aqf/`: shared core (data loading, quality rules, climatology, S2 features)
  - `Data Analysis/`: phase 1
  - `Forecasting/`: phase 2
- Git history: one commit per step.

**Phase 1: Data Analysis** (done, [`Data Analysis/REPORT.md`](<Data Analysis/REPORT.md>))
- Inventory and integrity: 438 five-day windows, one grid per sensor, SHA-256 provenance.
- Quality: **exact zeros are missing retrievals, not clean air**, so they're treated as missing.
- Trends: STL decomposition, seasonal Mann–Kendall / Sen slope, per-pixel trend maps with FDR correction.
- COVID lockdown effect, spatial hotspots, industrial hubs, and a Sentinel-2 summary.

**Phase 2: Forecasting** (complete: F0–F2, F4, F5 done, F3 skipped, [`Forecasting/PLAN.md`](Forecasting/PLAN.md), [`Forecasting/REPORT.md`](Forecasting/REPORT.md))
- Two tasks:
  - **A**: the next 5-day map.
  - **B**: the **150-day outlook** (30 leads).
- F0 forecastability: anomaly autocorrelation, cross-gas and neighbour predictive correlation, and Sentinel-2 land cover vs pollution with a spatial-block permutation test.
- F1 baselines: climatology, persistence and damped anomaly persistence. Metrics are skill, CRPS, interval coverage and block-bootstrap CIs.
- F2 learned models: LightGBM and Ridge per gas, pooled over pixels. Task A is a lead-1 model and task B a multi-horizon model with the lead as a feature. 27 features, including 30/60/90-day trailing anomalies. A no-leakage test covers the features, and there's a 2020 Q2 (lockdown) robustness check.
- F4: ablations (S2, neighbours, flags, trailing means, other gases, 2020 Q2), a **rolling-origin check over 2021, 2022 and 2023**, and an NE-monsoon gate. **The model choice was revised as a result.**
- F5: the frozen choice was scored **once on 2024** (fitted on 2019–2023), then everything was refitted on 2019–2024 for the **Jan–May 2025 outlook**: GeoTIFFs, CSVs, fan charts and hotspot series.

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

**But the F2 results don't hold outside 2023 (F4 rolling origin: train on the years before Y, score Y)**

| Year | NO2 5 d: Ridge / damped | NO2 10–30 d: Ridge / damped | NO2 35–150 d: Ridge / damped | CO 5 d: LightGBM / damped |
|---|---|---|---|---|
| 2021 | +2.0 / +3.9 | +0.9 / +0.6 | **−3.5** / −0.1 | −1.9 / **+18.8** |
| 2022 | −0.0 / +0.9 | +1.4 / −0.0 | +1.1 / −0.0 | +6.0 / +8.3 |
| 2023 | +12.7 / +6.6 | +6.7 / +1.2 | +5.6 / +0.5 | +15.5 / +9.4 |

- **The NO2 long-range skill depends on the year.** It is strong in 2023, small in 2022, and significantly *worse* than climatology in 2021.
- The learned 5-day models beat damped persistence only in 2023. Damped persistence is positive every year.
- Ablations (2023): the 30–90-day trailing means are the only input that matters for the long range. **Sentinel-2 adds nothing** (removing or shuffling it changes nothing). Neighbours and other gases add tenths of a percent.
- The NE-monsoon gate (climatology for Oct–Dec) doesn't help consistently, so it isn't adopted.

**Revised model choice (frozen for F5)**

| Gas | 5 days | 10–30 days | 35–150 days |
|---|---|---|---|
| NO2 | damped persistence | Ridge | climatology |
| CO | damped persistence | climatology | climatology |
| SO2 | climatology | climatology | climatology |

- The secondary candidates are NO2 Ridge at 5 days and at 35–150 days, and CO LightGBM at 5 days. They'll be scored once in F5 for the record.
- **Beyond about a month, the honest 150-day outlook is the seasonal climatology with calibrated 80 % bands** (coverage 0.76–0.83).
- Sentinel-2: mean NDVI vs long-term NO2 gives ρ = −0.53 (block-permutation p = 0.026). It explains *where* NO2 is high, not *when*.

### Final test: 2024, scored once with the frozen choice

| Gas | 5 d | 10–30 d | 35–150 d | 80 % coverage |
|---|---|---|---|---|
| NO2 | damped +3.1 % [−0.4, 6.1] | Ridge −0.5 % [−3.0, 3.0] | climatology (0) | 0.82 |
| CO | damped +3.8 % [−7.6, 10.1] | climatology (0) | climatology (0) | 0.85–0.87 |
| SO2 | climatology (0) | climatology (0) | climatology (0) | 0.81 |

- **2024 confirms the F4 revision and contradicts the 2023-only results.** The dropped NO2 Ridge long-range model is
  significantly worse than climatology in 2024 (−4.7 % [−7.3, −2.6]). The CO LightGBM 5-day model gives −1.4 %.
- The short-range gains are small: +3 to +4 %, with CIs touching 0. The intervals are calibrated.
- **Conclusion: beyond the seasonal cycle, predictability is small and year-dependent. The dependable 150-day forecast
  is the seasonal climatology with calibrated bands, plus a short-range correction.**

### Jan–May 2025 outlook (`Forecasting/results/F5_final/outlook_2025/`)
- 30 steps from 30 Dec 2024 to 29 May 2025, per gas: a 90-band GeoTIFF (mean, p10, p90; the model is named on every band),
  an area-mean CSV with 80 % bands, hotspot series (Manali NO2, Ennore SO2), fan charts and maps.
- Area-mean NO2 goes 35.5 µmol/m² (30 Dec) → about 40 (late Jan) → about 32 (late Mar low) → about 41 (late May), with 80 % bands of ±12–13.

## 3. Not done yet
- Nothing from the plan remains. Possible extensions:
  - add meteorology (wind, boundary-layer height, rain), the most likely source of extra skill
  - re-score when 2025 data arrive
  - write up the final report and slides

## 4. Risks and limitations
- **Six years is short.** Trend tests have low power, and only one validation year and one test year exist.
- **SO2 is noisy.** 39 % of its pixels are zeros, 21 % of its windows are unusable, and nothing beats climatology. Its outlook will be climatology only.
- **2023 (validation) was the highest-NO2 year** and flattered the learned models. The rolling-origin check exposed this, and the choice was revised.
- **Model selection used 2021–2023**, so the single 2024 score is the only unbiased estimate. It agrees with the revised choice.
- **NO2 is weak in the NE monsoon (Oct–Dec)**, and a simple gate didn't fix it.
- **Validation numbers drive model choice** and are therefore optimistic. Only the single 2024 score is unbiased.
- **The old pipeline exists only in the system trash.** Don't empty the trash if any of it might still be needed.

## 5. Commands (run from `Scratch/`)

```bash
uv sync
uv run pytest -q
uv run aqf-analysis all      # phase 1, about 1 min, results in "Data Analysis/results/"
uv run aqf-forecast all      # phase 2 F0 + F1 + F2 + F4, about 4 min, results in Forecasting/results/
uv run aqf-forecast outlook  # refit on 2019-2024 and write the Jan-May 2025 outlook (the 2024 test is locked)
```
