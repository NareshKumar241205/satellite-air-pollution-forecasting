# Forecasting: findings so far (F0 forecastability, F1 baselines)

To reproduce, run `uv run aqf-forecast all` (about 40 s). Results go to `results/F0_forecastability/` and `results/F1_baselines/`.
The plan and the remaining phases are in [`PLAN.md`](PLAN.md). Every number below comes from fits on the
**train years 2019–2022** and is scored on the **validation year 2023**. The 2024 test year hasn't been touched.

- Task A is the next 5-day map (lead 1).
- Task B is the 150-day outlook (leads 1–30).
- "Skill" is the MSE skill against the per-pixel seasonal climatology: 1 − MSE_model / MSE_climatology. Positive means
  better than the seasonal average. The 95 % CIs come from a moving-block bootstrap over forecast origins.

## F0: How forecastable are the anomalies?

Anomaly means the departure from the per-pixel seasonal climatology fitted on train.

| | NO2 | CO | SO2 |
|---|---|---|---|
| Autocorrelation at +5 days, single pixel | 0.13 | **0.35** | −0.01 |
| Autocorrelation at +5 days, area mean | 0.38 | **0.45** | 0.21 |
| Area-mean autocorrelation significant up to | 30 days | 10 days | 5 days |
| Pixel autocorrelation drops below 0.1 from | +10 days | +15 days | +5 days |

(`acf.csv`, `acf_reach.csv`, `fig_acf.png`)

- **The memory of the anomalies is short.** At pixel level it is gone after 1–3 windows, so day-to-day weather and
  retrieval noise dominate a single 5 km pixel. Area means keep more memory because averaging removes the noise.
- **NO2 has a slow component.** The area-mean NO2 anomaly stays positively correlated (+0.1 to +0.23) out to about
  90 days, at or above the noise band. Whole periods run high or low, 2020 for example. The baselines below can't use this; F2 will
  test trailing 30/90-day mean anomalies as inputs. Caveat: the 2020 lockdown low sits in the train years and may
  produce part of this pattern.
- **Other gases add almost nothing at the pixel level.** The partial correlation of gas A at t with gas B at t+5 days,
  given B's own value, is at most 0.026 (`cross_gas_partial.csv`). The monthly NO2–SO2 co-variation found in the analysis
  doesn't carry over into predictive information at 5-day pixel scale.
- **The neighbourhood adds a little.** The 3×3 neighbour mean predicts the next window beyond the pixel's own value with
  partial r = 0.09 (NO2) and 0.14 (CO) at +5 days, fading by +15 days. SO2 is ≈ 0 (`neighbour_partial.csv`).
- **Sentinel-2 land cover vs the long-term mean**, over the 92 mostly-land S5P pixels, with a spatial-block permutation test (`s2_vs_longterm_mean.csv`):

| S2 feature (train years) | NO2 ρ | p (block perm.) | p (naive) |
|---|---|---|---|
| built-up share (NDBI > 0) | 0.34 | 0.14 | 0.0009 |
| mean NDBI | 0.38 | 0.09 | 0.0002 |
| mean NDVI | **−0.53** | **0.026** | <0.0001 |

  Greener pixels have less NO2. That survives the spatial permutation test. Built-up share points the same way
  but isn't significant once spatial autocorrelation is respected. The naive p-values overstate the evidence by
  about 100×. CO and SO2 show no relation. **S2 can at most help explain *where* NO2 is high, which the per-pixel climatology
  already captures, not *when*.** So any S2 benefit in forecasting is expected to be small. F4 tests this directly (S2 off and S2 shuffled).

## F1: Baselines

### Task A, the next 5-day map (`taskA_lead1_val.csv`)

| Gas | Model | Pixel skill [95 % CI] | Area-mean skill | 80 % interval coverage |
|---|---|---|---|---|
| NO2 | persistence | −39 % [−58, −14] | −1 % | 0.78 |
| NO2 | **damped anomaly persistence** | **+6.6 % [+4.3, +9.9]** | **+11 %** | 0.79 |
| CO | persistence | −43 % [−73, −12] | −21 % | 0.78 |
| CO | **damped anomaly persistence** | **+9.4 % [+0.1, +21]** | **+16 %** | 0.79 |
| SO2 | persistence | −83 % [−102, −68] | −40 % | 0.78 |
| SO2 | damped anomaly persistence | −0.3 % [−0.5, −0.1] | −0.3 % | 0.79 |

Lead-1 RMSE: NO2 19.6 vs 20.2 µmol/m² for climatology, CO 3.71 vs 3.89 mmol/m², SO2 175 vs 175 µmol/m².

- **Persistence (repeating the last map) is much worse than the seasonal average**, because a single 5-day map is noisy.
- **Damped anomaly persistence is the baseline to beat.** It starts from climatology and adds the last anomaly, shrunk by the
  measured autocorrelation. It beats climatology significantly for NO2 and CO. For SO2 nothing beats climatology.

### Task B, the 150-day outlook (`taskB_summary_val.csv`, `fig_baselines_by_lead_val.png`)

| Gas | Damped persistence beats climatology (unbroken CI > 0) for | Mean skill over leads 1–30 |
|---|---|---|
| NO2 | 6 leads = **30 days** | +0.8 % |
| CO | 1 lead = 5 days | +0.2 % |
| SO2 | none | −0.1 % |

- **Beyond about a month, the best available forecast is the seasonal climatology.** For leads 35–150 days, every
  baseline's skill is indistinguishable from 0. For CO, the climatology alone is already a good long-term outlook, because
  its seasonal cycle explains 83 % of the variance.
- **The intervals are calibrated.** The 80 % intervals cover 78–80 % of validation observations at every lead for climatology
  and damped persistence. They are built from out-of-sample residuals (leave one train year out).
  Persistence intervals are miscalibrated (NO2 71–78 %, SO2 76–78 %, CO up to 87 %), and its CRPS grows with lead.

## Implications for F2 and later

1. Task A has room to improve on +6.6 % (NO2) and +9.4 % (CO). The candidate inputs, in order of F0 evidence, are:
   - the 3×3 neighbourhood
   - **trailing multi-window mean anomalies (30–90 days)**
   - season interactions
   - other gases and S2 last, since their F0 signal is weak
2. For task B, a learned model must beat climatology beyond 30 days to add anything. The slow NO2 component is the
   only lead. Expect the outlook to be **climatology plus a decaying anomaly correction**, and report it that way.
3. **SO2 isn't forecastable beyond its seasonal climatology** with these data. Report its outlook as climatology with
   calibrated intervals.
4. **F3 (ConvLSTM) isn't justified by F0 so far.** The spatial and temporal dependence it would model is weak
   (neighbour partial r ≤ 0.14, pixel memory under 15 days). It stays conditional on F2.

---

# F2: Learned models (LightGBM + Ridge)

To reproduce, run `uv run aqf-forecast learned` (about 3 min). Results are in `results/F2_learned/`.

**Setup.**
- There is one model per gas, pooled over the 144 pixels:
  - task A: a dedicated lead-1 model
  - task B: one direct multi-horizon model with the lead as a feature
- Features (27): own-pixel anomaly lags 0–5 with missing flags, 3×3 neighbour lags 0–2, area-mean anomaly,
  **trailing 30/60/90-day area and pixel mean anomalies**, target-window climatology and season, pixel position, and 5 static S2
  features. Other gases are left out (F0: partial r ≤ 0.026).
- Training uses train-year targets only. LightGBM early stopping and the Ridge α use validation 2023, so **validation
  scores are optimistic**.
- Intervals come from leave-one-train-year-out residuals, the same method as the baselines.
- `tests/test_forecast.py::test_no_lookahead_in_f2_features` checks that no feature changes when the future is altered.

## Task A: next 5-day map (`A_taskA_lead1_val.csv`)

| Gas | Best model | Pixel skill vs climatology [95 % CI] | vs damped persistence [95 % CI] | Area-mean skill |
|---|---|---|---|---|
| NO2 | **Ridge** | **+12.7 % [3.6, 22.9]** | +6.5 % [−1.3, 15.0] | **+35 %** |
| NO2 | LightGBM | +6.1 % [0.0, 14.1] | −0.6 % [−5.3, 5.7] | +16 % |
| CO | **LightGBM** | **+15.5 % [2.9, 29.3]** | +6.8 % [−1.6, 14.4] | **+23 %** |
| CO | Ridge | +12.1 % [1.0, 26.3] | +3.0 % [−0.3, 6.7] | +18 % |
| SO2 | LightGBM | +0.2 % [−0.1, 0.6] | +0.4 % [0.2, 0.8] | +1 % |

For reference, damped persistence scored +6.6 % for NO2 and +9.4 % for CO.

- **The learned models roughly double the skill of the best baseline for NO2 and CO**, and lift the area-mean skill to 23–35 %.
- The improvement over damped persistence is **not yet statistically significant**: both CIs include 0.
  A single validation year limits the power. F5 (2024) is the real check.
- SO2 is still not forecastable. The tiny gain is significant but practically irrelevant.

## Task B: 150-day outlook (`B_taskB_summary_val.csv`, `fig_f2_by_lead_val.png`)

| Gas | Model | Beats climatology (unbroken 95 % CI > 0) | Beats damped persistence (CI > 0) | Mean skill, leads 1–30 | Area-mean skill at 5 / 30 / 60 / 90 / 120 / 150 d |
|---|---|---|---|---|---|
| NO2 | **Ridge** | **leads 1–26 = 130 days** | **leads 2–26** | **+5.9 %** | 24 / 19 / 23 / 22 / 13 / 10 % |
| NO2 | LightGBM | all 30 leads = 150 days | most of leads 6–30 | +3.2 % | 6 / 7 / 7 / 6 / 10 / 10 % |
| CO | Ridge | 10 days | leads 6–11 | +1.2 % | 7 / 6 / 5 / 2 / −1 / −4 % |
| CO | LightGBM | none | none | +0.7 % | ≈ 0–2 % |
| SO2 | all | none | none | ≈ 0 | ≈ 0 |

- **The 150-day outlook now carries real information for NO2.** The learned models beat the seasonal climatology for
  4–5 months, where the baselines had nothing after 30 days. The source is the slow, area-wide NO2 component that F0 found.
  The top LightGBM features after season are the **60- and 90-day trailing area-mean anomalies** (`feature_importance.csv`).
- For CO, the seasonal climatology remains the outlook beyond about 10 days. For SO2 it is the outlook at every lead.

### Does the NO2 long-lead gain survive without the 2020 lockdown? (`check_2020q2.csv`)

| NO2, mean skill at leads 35–150 days | All train years | Train without Apr–Jun 2020 |
|---|---|---|
| Ridge | +5.6 % | **+6.9 %** |
| LightGBM | +3.3 % | +3.5 % |
| Damped persistence | +0.5 % | — |

**Yes: the gain survives, and slightly grows,** when the lockdown months are dropped from training. The lockdown
does inflate the raw area-mean autocorrelation (30-day ACF 0.14 → 0.07 without it), but the models' skill doesn't
depend on it.

### Caveat
2023 was the highest-NO2 year in the record (+14 % vs 2019), so "the recent anomaly persists" was a good bet that year.
2024 was 15 % lower than 2023. The single 2024 test (F5) will show whether the long-lead skill holds in a year that
changes direction.

## Seasons, intervals and features
- **By season** (`season_breakdown_val.csv`): the learned models gain most in the pre-monsoon and SW monsoon seasons
  (NO2 task A: Ridge +25 % and +16 %; CO: LightGBM +20 % and +33 %). They are **worse than climatology in the NE monsoon
  (Oct–Dec)** for NO2, by 1–9 %. That's a known weak spot.
- **Intervals:**
  - Residual-based 80 % intervals cover 0.76–0.83 (mostly 0.78–0.80), so they're close to calibrated.
  - Quantile-LightGBM bands under-cover at 0.73–0.77, so **F5 will use residual intervals**.
  - CRPS is lowest for Ridge (NO2) at almost every lead.
- **Features:** season and the trailing area means dominate. The neighbourhood contributes a little to task A. The static S2
  features aren't among the top features for any gas, consistent with F0. F4 tests this directly.
- **LightGBM early-stopped after 2–33 rounds.** The signal is weak and close to linear, which is why Ridge matches or beats it.

## Decision
- **F3 (ConvLSTM) isn't justified.** The useful signal is area-wide and slow, near-linear, and not fine-scale
  spatio-temporal. Pixel and neighbour features add little.
- **Candidate models frozen for F5:**

  | Gas | Task A | Task B |
  |---|---|---|
  | NO2 | Ridge | Ridge |
  | CO | LightGBM | climatology beyond 10 days |
  | SO2 | climatology | climatology |

  All use residual intervals. F4 ablations come first, to confirm which inputs matter.

---

# F4: Ablations, rolling-origin robustness, seasonal gate

To reproduce, run `uv run aqf-forecast ablations` (about 30 s; it needs F2 first). Results are in `results/F4_ablations/`.
Hyperparameters are frozen from F2: Ridge α and LightGBM rounds. Brackets are 95 % block-bootstrap CIs.

## 1. Ablations on validation 2023 (`ablations_val.csv`, `fig_ablations.png`)

One change at a time. The table shows skill relative to the full model; negative means the change hurts.

| Change | NO2 task A (5 d) | NO2 task B, 35–150 d | CO task A (5 d) |
|---|---|---|---|
| no Sentinel-2 features | +0.0 [−0.0, +0.1] | −0.0 | −1.8 [−5.7, +1.7] |
| Sentinel-2 shuffled across pixels | +0.0 | −0.0 | +0.0 |
| no 3×3 neighbours | **−0.4 [−0.6, −0.1]** | +0.0 | −0.1 |
| no missing-data flags | +0.3 | **−0.4 [−0.7, −0.2]** | −0.3 |
| **no 30/60/90-day trailing means** | −2.4 [−5.7, +0.9] | **−4.3 [−6.2, −2.1]** | +0.3 |
| add the other gases | −0.2 | +0.5 | +0.5 |
| train without Apr–Jun 2020 | −0.5 | **+1.3 [+0.9, +1.7]** | +1.4 |

- **The trailing (30–90-day) mean anomalies are the only input that matters for the long range.** Without them the NO2
  35–150-day skill falls from +5.6 % to +1.6 %.
- **Sentinel-2 adds nothing to forecasting.** Removing or shuffling it doesn't change skill. S2 describes *where* NO2
  is high, which the per-pixel climatology already holds (F0).
- Neighbours, missing-data flags and other gases are worth at most a few tenths of a percent.

## 2. Rolling origin: does the F2 choice hold in other years? (`rolling_origin.csv`, `fig_rolling_origin.png`)

For each year Y, everything is refit on the years before Y and Y is scored. 2021 trains on two years and 2022 on three, so the
learned models get less data in those years. Skill is against the seasonal climatology.

| Year | NO2 task A: Ridge / damped | NO2 10–30 d: Ridge / damped | **NO2 35–150 d: Ridge / damped** | CO task A: LightGBM / damped |
|---|---|---|---|---|
| 2021 | +2.0 / +3.9 | +0.9 / +0.6 | **−3.5 [−6.7, −0.9]** / −0.1 | −1.9 / **+18.8** |
| 2022 | −0.0 / +0.9 | **+1.4 [+0.3, +2.3]** / −0.0 | +1.1 [+0.3, +1.9] / −0.0 | +6.0 / +8.3 |
| 2023 | **+12.7** / +6.6 | **+6.7** / +1.2 | **+5.6** / +0.5 | **+15.5** / +9.4 |

**The F2 conclusions don't hold outside 2023. Plainly:**
- **The NO2 long-range (35–150 d) skill isn't robust.** It is +5.6 % in 2023 and only +1.1 % in 2022. In 2021 it is
  **significantly worse than the seasonal average (−3.5 %)**. The slow NO2 component exists, but whether it helps depends on the year.
- **The task A learned models don't beat damped persistence reliably.**
  - NO2: Ridge beats damped persistence only in 2023.
  - CO: LightGBM loses to damped persistence in 2021 and 2022. In 2021 it is 25 % worse (significant).
- **What does hold every year:**
  - Damped persistence is positive at 5 days in every year for both gases.
  - NO2 Ridge at 10–30 days is ≥ damped persistence in all three years (significant in 2022 and 2023).

## 3. NE-monsoon gate (`seasonal_gate.csv`)

The gate replaces the model forecast with climatology for targets in Oct–Dec.
- It helps in some year/gas combinations (CO task A: +3.6 in 2021, +0.8 in 2022) and hurts in others
  (CO 2023 −1.5; NO2 2021 and 2022: −0.5 to −1.4).
- **It doesn't help consistently, so it isn't adopted.** The Oct–Dec weakness is real, but a hard switch to climatology doesn't fix it.

## Revised model choice (frozen for F5)

| Gas | Task A (5 d) | Task B, 10–30 d | Task B, 35–150 d |
|---|---|---|---|
| NO2 | **damped persistence** (positive every year) | **Ridge** (≥ damped in all 3 years) | **climatology** |
| CO | **damped persistence** | climatology | climatology |
| SO2 | climatology | climatology | climatology |

- All of these use the out-of-sample residual intervals, which are close to calibrated.
- Two models are kept as **secondary candidates** and scored once in F5 next to the primary choice, so the report can say how they did. They don't change the primary forecast:
  - NO2 Ridge at 5 days and at 35–150 days
  - CO LightGBM at 5 days
- **Bottom line for the 150-day outlook.** Beyond about a month, the honest forecast is the seasonal climatology with calibrated
  bands. The satellite record has information about *where* and *in which season* pollution is high, but not about
  whether the next 2–5 months will be above or below normal.

---

# F5: The single 2024 test and the Jan–May 2025 outlook

To reproduce, run `uv run aqf-forecast final` (about 1 min). Results are in `results/F5_final/`.

**How the test was run**
- The frozen choice from F4 was scored **once** on 2024.
- Every component was fitted on 2019–2023 (climatology too) with the hyperparameters frozen in F2.
- Nothing was changed after seeing these numbers. A lock file (`TEST_SCORED.json`) blocks a silent second scoring.
- One re-run was made with `--rescore-test`, only to fix a table bug: lead-1-only models were listed at leads they don't
  forecast. The numbers were identical (max difference 0.0).

## 2024 test results (`test_summary.csv`, `fig_test_2024_by_lead.png`)

Skill is against the 2019–2023 seasonal climatology. Brackets are 95 % block-bootstrap CIs.

| Gas | Leads | Frozen model | 2024 skill | 80 % interval coverage | Alternatives (for the record) |
|---|---|---|---|---|---|
| NO2 | 5 d | damped persistence | +3.1 % [−0.4, +6.1] | 0.82 | Ridge (A) +3.9 % [−0.4, +9.1] |
| NO2 | 10–30 d | Ridge (B) | −0.5 % [−3.0, +3.0] | 0.82 | damped +0.0 % |
| NO2 | 35–150 d | climatology | 0 (reference) | 0.82 | **Ridge (B) −4.7 % [−7.3, −2.6]** |
| CO | 5 d | damped persistence | +3.8 % [−7.6, +10.1] | 0.85 | LightGBM (A) −1.4 % [−24, +13] |
| CO | 10–150 d | climatology | 0 (reference) | 0.87 | damped −1.2 % to +0.1 % |
| SO2 | 5–150 d | climatology | 0 (reference) | 0.81 | damped −0.1 to −0.4 % |

**What 2024 says, plainly:**
- **2024 contradicts the 2023 validation results for the learned models, and confirms the F4 revision.**
  - The NO2 Ridge long-range model, which scored +5.6 % on 2023, is **significantly worse than climatology in 2024 (−4.7 %)**.
  - The CO LightGBM 5-day model, which scored +15.5 % on 2023, gives −1.4 %.
  - Keeping them would have made the outlook worse.
- **The short-range gains are small in 2024.**
  - Damped persistence at 5 days gives +3.1 % (NO2) and +3.8 % (CO). Both are positive, but the CIs include 0.
  - The NO2 Ridge at 10–30 days, the one learned component kept, has no skill in 2024 (−0.5 %).
- **The intervals hold.** 80 % bands cover 81–87 % of 2024 observations, so they're calibrated to slightly conservative (CO).
- **Overall:** across 5–150 days, the frozen choice performs the same as the seasonal climatology (NO2 +0.0 %, CO +0.1 %,
  SO2 0). Predictability beyond the seasonal cycle is small and year-dependent. The **seasonal climatology with calibrated
  uncertainty is the dependable forecast** for this record, and the main deliverable is that climatology plus a short-range
  correction, with honest bands.

## Jan–May 2025 outlook (`outlook_2025/`)

Everything was refitted on 2019–2024, with the climatology recomputed on those six years. The forecast is issued from the
25 Dec 2024 window for 30 steps, **30 Dec 2024 to 29 May 2025**. The model per lead is the frozen choice:
- NO2: damped at 5 d, Ridge at 10–30 d, climatology beyond
- CO: damped at 5 d, climatology beyond
- SO2: climatology

| Gas | 30 Dec–4 Jan | late Jan | late Mar (seasonal low) | late May | Unit |
|---|---|---|---|---|---|
| NO2, area mean | 35.5 [24.5, 48.0] (damped) | 40.5 [28.5, 53.0] (Ridge, +30 d) | 32.2 [20.0, 45.6] | 41.5 [29.1, 54.8] | µmol/m² |
| CO, area mean | 41.0 [37.5, 45.0] (damped) | 39.9 [36.1, 44.5] | 43.5 [39.7, 48.1] (seasonal high) | 36.0 [32.3, 40.6] | mmol/m² |
| SO2, area mean | 197.8 [124.3, 272.2] | 217.9 [144.8, 293.8] | 157.1 [83.4, 236.5] | 194.6 [119.9, 275.3] | µmol/m² |

Values are the mean with the 80 % band in brackets.

**Files**
- `outlook_<GAS>_2025.tif`: 90 bands each. Bands 1–30 are the mean, 31–60 the 10 % quantile, 61–90 the 90 % quantile, in mol/m².
  Each band's description and tags name the lead, the window and the **model used**.
- `outlook_area_mean.csv`: one row per gas and lead, with the model column.
- `outlook_hotspots.csv`: Manali/CPCL NO2 and Ennore SO2 hub pixels.
- `fig_outlook_fan_area.png`, `fig_outlook_hotspots.png`, `fig_outlook_maps.png`, plus `README.md`.

At Manali, NO2 is forecast to peak at about 84 µmol/m² in late January and fall to about 38 in early April.
That is the hotspot's typical seasonal pattern, with an 80 % band of roughly ±20 µmol/m².
