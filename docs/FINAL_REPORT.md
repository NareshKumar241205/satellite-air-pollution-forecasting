# Final report: land-use change and Sentinel-5P NO2 over Tamil Nadu, with a 5-day to 150-day forecast

Branch `adr-001-phase1`, 12 Oct 2026. Every number below comes from files in `Land Change/results/`. Decision rules were committed before the results they govern:
- `Land Change/SITE_SELECTION.md`
- `Land Change/DIAGNOSTICS_PLAN.md`
- `Land Change/DECISION_RULES.md`

## 1. Question and answer

*Other gases (CO, SO2) are summarised in answer 4 below and detailed in §4b.*

**Question.** Does land-use change (built-up growth) measurably change tropospheric NO2 in Sentinel-5P? Can a
satellite NO2 forecast for Tamil Nadu cities use land change, and how far ahead is it skilful?

**Answer.**
1. **No land→NO2 effect is detectable.** On the best available design, the smallest effect that could have been detected (80 % power) is **0.73 % NO2 per percentage point of built-up gain**. Every estimate is near zero; for example, own built-up gives −0.004 %/pp [−0.08, 0.09]. The design checks fail, though: the lead placebo is non-zero, and it is driven by the 2020 lockdown year. So the null can't be read as "no effect" either. The thesis result is the detection limit.
2. **The land term in the forecast is therefore zero (B = 0),** under a rule fixed in advance.
3. **Forecasts beat climatology only for the next 5-day window.**
   - Validation 2023–24: **+9.4 %** MSE skill [6.4, 13.4] (Ridge).
   - Test 2025, scored with frozen settings: **+3.6 %** [−1.3, 6.8].
   - From 10 days to 150 days, skill is ≈ 0, so the dependable long-range forecast is the per-pixel seasonal climatology with calibrated bands. The 10–90 % band covered 81 % of 2025 pixel values.
4. **Other gases.** CO is forecastable for the next window (**+25.8 %** [12.0, 38.6] on the 2025 test, but bands too narrow: 0.74 coverage); SO2 isn't (≈ 0). Land has no detectable effect on CO (below ≈ 0.02 %/pp, a precise null) and SO2's land test is uninformative (limit 4.7 %/pp).

## 2. Data

| Item | Source | Notes |
|---|---|---|
| NO2 | `COPERNICUS/S5P/OFFL/L3_NO2` tropospheric column, 5-day masked means on 0.01° (TNAQ export, `D:\CLG\SEM-7\FINAL PROJECT\Dataset\Raw`, read-only), 2018-07 → 2024-12 | 2025: 73 windows downloaded with the same method into `data/tnaq_extra_v2` (§7) |
| Analysis grid | 4×4 blocks = 0.04° (≈ 4.4 km), 12×12 blocks per site, **18 sites**, 2,568 blocks | Boxes don't overlap (a test enforces it) |
| Land | Dynamic World V1 same-season (Jan–Apr) mode label, 10 m → built-up fraction per block and year, 2019–2025 | Alternatives: DW `built` probability (trailing 3-year median), Open Buildings Temporal |
| Weather | ERA5 hourly (`ECMWF/ERA5/HOURLY` in Earth Engine): 06–09 UTC, site box + 0.25°, per 5-day window, 2019–2025 | Same fields as the CDS product; no CDS account needed |
| Processor versions | Earth Engine `PROCESSOR_VERSION` per image | The Jan–Apr seasons use 1.2.2/1.3.x, 1.3.2, 1.4.0, 2.3.1, 2.4.0/2.5.0, 2.6.0 and 2.8.0 for 2019–2025; v1→v2 falls between 2021 and 2022 |

**Splits.**
- Land analysis: 2019–2023 only.
- Forecast: train 2019–2022, validation 2023–2024, test 2025.

**Out of reach.** The old Chennai `Dataset/data` and the F0–F5 run outputs live on the Linux machine, not on this PC. Their 2024 test record is in `Forecasting/results/F5_final/` and is unchanged.

## 3. Land change → NO2

### 3.1 Design and gate G1

- **Unit:** one block × one year. y = 100·log(Jan–Apr mean NO2), L = built-up %.
- **Gate G1:** the pooled long difference 2019→2023 must be able to detect 0.5 %/pp at 80 % power, after removing land-label noise and spatial autocorrelation.
- **Site selection:** sites could be added only by DW built-up growth. The rule was frozen before any slope. 0 of 10 candidate tiles reached the 1.5 pp threshold.

| Built-up layer | Detectable effect (%/pp) | G1 at 0.5 | at 1.0 |
|---|---|---|---|
| DW mode label | **0.73** | fail | pass |
| DW probability, trailing 3-y median | 1.52 | fail | fail |
| Open Buildings Temporal | 4.17 | fail | fail |
| 1 km S2 index proxy (early check) | ∞: it measured greenness, r = −0.99 with vegetation change | fail | fail |

The smoother layers move less from year to year, so they carry less signal.

### 3.2 Exploratory panel (C2)

- **Model:** pixel + site×year fixed effects, own and neighbour built-up.
- **Inference:** a 4×4-block spatial bootstrap (500 draws), plus a spatial placebo, a lead/lag test and the long difference.
- **Outcomes:** raw and ERA5-deweathered NO2.

| Layer / outcome | own L | long difference | lead L(t+1), should be 0 | neighbour WL |
|---|---|---|---|---|
| DW label, raw | −0.004 [−0.08, 0.09] | −0.03 [−0.15, 0.10] | **+0.31 [0.14, 0.45]** | +0.27 [0.05, 0.48] |
| DW label, deweathered | −0.002 [−0.08, 0.07] | +0.00 [−0.13, 0.14] | **+0.35 [0.19, 0.46]** | +0.41 [0.20, 0.60] |
| DW prob. 3-y, raw | −0.03 [−0.16, 0.11] | −0.07 [−0.28, 0.13] | +0.24 [−0.27, 0.77] | +0.17 [−0.22, 0.58] |
| DW prob. 3-y, deweathered | −0.02 [−0.16, 0.11] | +0.00 [−0.22, 0.19] | **+0.54 [0.02, 1.02]** | +0.69 [0.27, 1.07] |
| Open Buildings, raw | −0.39 [−1.26, 0.87] | +0.13 [−1.05, 1.47] | **+2.75 [0.36, 4.87]** | +1.66 [−0.56, 4.40] |

### 3.3 Why the lead test fails (`Land Change/results/D_lead_diagnostics/`)

**Found:**
- **The lockdown year.** Dropping outcome year 2020 removes the lead: 0.31 → **0.01** (−97 %). Jan–Apr 2020 contains the COVID lockdown (from 25 March), whose NO2 drop varied across space. The likeliest reading is that pixels which later gained built-up area had different lockdown drops. D4 shows only that 2020 carries the lead, not the mechanism.


**Ruled out:**
- **Aerosol, brightness, cloud and scene counts (D1, all of 2019–2024).** Within site and year, aerosol optical depth tracks the built-up label strongly (r = 0.56) and the label lead (0.59), as does scene brightness (0.51/0.50). But the correlation with the NO2 residual is ≈ 0 (AOD +0.02, −0.02 next year; brightness +0.02), so they don't meet the pre-set flag (|r| ≥ 0.1 with both the land lead and NO2). An earlier reading from 2019–2021 only gave −0.12 for next-year AOD; with the full six years it is −0.02, and that statement is withdrawn. Cloud probability and the number of Dynamic World images show no relation to either (|r| < 0.1).
- **Single sites:** leave-one-site-out moves the lead by at most ±24 %.
- **Box edges:** dropping the outer ring changes it by −11 %.
- **Shared spatial structure:** the neighbour cross-correlation of the residuals is 0.02.

**Partly explained:**
- **Pixel trends:** pixel-specific trends shrink the lead to 0.11.
- **Level dependence:** yearly label change is larger in high-NO2 pixels in the first three pairs: +0.8, +1.2 and +2.0 pp per 100 log-units for 2019→20, 2020→21 and 2021→22. It is −0.7 for 2022→23. That is consistent with growth on the polluted urban fringe, but it isn't stable across years.

**Post-hoc D7: the panel without the 2020 outcomes** (declared in `DIAGNOSTICS_PLAN.md` before it was run; it can't change rule B). Outcome years 2021–2022 remain for the lead test.

| Layer | own L, %/pp [95 % CI] | lead L(t+1) | neighbour WL |
|---|---|---|---|
| DW label, raw | +0.001 [−0.08, 0.10] | +0.009 [−0.27, 0.25] | +0.25 [0.02, 0.49] |
| DW label, deweathered | +0.011 [−0.07, 0.09] | +0.010 [−0.21, 0.24] | +0.35 [0.14, 0.54] |
| DW prob. 3-y, raw | −0.04 [−0.19, 0.12] | +0.25 [−0.49, 1.03] | +0.23 [−0.26, 0.70] |
| Open Buildings, raw | −0.77 [−1.78, 0.61] | **−3.56 [−6.49, −0.43]** | +0.84 [−1.47, 4.02] |

Without 2020 the lead is zero on the main layer, which confirms the lockdown reading. The neighbour term stays positive on the DW-label layer (spillover or residual confounding; this design can't separate them), and Open Buildings stays unusable.

**Rule B (pre-registered).**
- Steadier layer: DW probability, the only one clean on the primary (raw) outcome.
- G1 fails at 0.5 %/pp, so **B = 0**.

## 4. Forecasting (C3–E)

**Target.** The block-mean NO2 of the next window, as an anomaly from a per-pixel 3-harmonic climatology.

**Features.** All known at the origin:
- own lags 0–5 and missing flags;
- the 3×3 neighbour mean and the site mean, lags 0–2;
- trailing means over 6, 12 and 18 windows;
- day-of-year harmonics;
- ERA5 scalars of windows t and t−1;
- the built-up fraction of the last complete season.

### 4.1 Lead 1 (next 5 days), validation 2023–2024

| Model | MSE skill vs climatology | RMSE (µmol/m²) | MAE | pattern r |
|---|---|---|---|---|
| climatology | 0 | 10.20 | 6.70 | – |
| damped persistence (φ₁ = 0.15) | +5.4 % [3.7, 7.1] | 9.93 | 6.50 | 0.13 |
| **Ridge (α = 1000), chosen** | **+9.4 % [6.4, 13.4]** | 9.71 | 6.35 | 0.11 |
| LightGBM (48 rounds) | +6.5 % [4.0, 9.1] | 9.87 | 6.44 | 0.10 |
| small CNN | not run: the PyTorch download was stopped | | | |

Ridge validation skill by IMD season: winter +18 %, pre-monsoon +10 %, SW monsoon +10 %, NE monsoon +5 %.

The chosen α = 1000 is the largest value in the frozen grid {0.1 … 1000}. Stronger shrinkage might do slightly better, but the grid wasn't widened after seeing it, so it stays.

**Ablations** (LightGBM, validation skill):

| Variant | Skill |
|---|---|
| full | +6.49 % |
| ERA5 off | +5.81 % |
| land off | +6.48 % |
| land shuffled within site | +6.49 % |

Weather helps a little; **land adds nothing**.

### 4.2 Long range and fusion

- **Fusion:** F = climatology + B + g(h)·Â(h), with g(h) fitted per lead on validation and monotone in [0, 1]. Here g ≈ 1 up to 55 days and 0.85 from 60 days.
- **Skill by lead:** damped persistence falls to +1.7 % at 10 days, +1.0 % at 15 days, and below +1 % after that.
- **Long-range forecast:** beyond the first window, the forecast is effectively climatology. See `Land Change/results/figures/fig_skill_vs_lead.png`.

### 4.3 Test 2025, frozen settings

| Lead | Fused skill | RMSE (clim RMSE) | 10–90 % coverage |
|---|---|---|---|
| 5 days | **+3.6 % [−1.3, 6.8]** | 9.31 (9.48) | 0.81 |
| 10 days | +0.05 % | 9.47 | 0.81 |
| 30 days | +0.5 % | 9.46 | 0.81 |
| 60–150 days | ≈ 0 (mean over leads 2–30: −0.05 %) | 9.47–9.48 | 0.81 |

- **Pattern correlation** at lead 1 is 0.09.
- **Site-mean fan chart** from the last 2024 window: 78 % of observed site means fall inside the 10–90 % band.
- **Bias:** Jan–May 2025 site means ran a median 7 % below the 2019–2024 climatology.

## 4b. The other gases: CO and SO2 (rule G, committed before any result)

Same sites, blocks, splits, features and models. Land weight stays 0 for both gases. Each gas was scored on 2025 exactly once, using its own 2025 data (73 windows per gas, reproduced bit for bit on two checked windows each). Units: mmol/m² for CO and µmol/m² for SO2. Results are under `Land Change/results/gas/<gas>/`, and the comparison figure is `fig_skill_by_gas_2025.png`.

### Forecast skill

| | NO2 | CO | SO2 |
|---|---|---|---|
| Validation 2023–24, lead 1, best model | Ridge +9.4 % [6.4, 13.4] | Ridge **+16.3 % [5.2, 27.0]** | none: all models −0.7 % to 0.0 % |
| (damped persistence / LightGBM) | +5.4 % / +6.5 % | +13.2 % / +13.1 % | −0.1 % / 0.0 % |
| **2025 test, lead 1** | **+3.6 % [−1.3, 6.8]** | **+25.8 % [12.0, 38.6]** | +0.1 % [0.0, 0.2] |
| 2025 test, lead 2 (10 days) | +0.05 % | +5.0 % | −0.1 % |
| 2025 test, leads ≥ 6 (30 days+) | ≈ 0 | ≈ 0 (−1.8 % to +0.7 %) | ≈ 0 |
| 10–90 % band coverage, 2025, lead 1 | 0.81 | 0.74 (0.68–0.70 beyond) | 0.94 |
| Reading (rule G) | skill in week one only | **real skill in the next window**; climatology after | **no skill beyond climatology** |

- **CO.** Its column is dominated by slow, large-scale variability, so recent anomalies persist: damped persistence alone gives +13 % at lead 1. 2025 was a more variable year than validation (climatology RMSE 3.66 against 3.11 mmol/m²), which raised the test skill but left the bands too narrow (coverage 0.74 against 0.80 nominal).
- **SO2.** The retrievals are very noisy (RMSE ≈ 200 µmol/m² around the climatology, from a wide negative tail allowed by the −0.005 mol/m² floor), and the signal is mostly point sources. Validation skill is −0.7 % to 0.0 %, so by the reading rule SO2's forecast is the climatology. The 0.94 coverage reflects bands that are too wide for heavy tails.
- **Ablations (LightGBM, validation).** Weather helps CO (+13.1 % with it, +11.8 % without); neither SO2 nor CO uses land. The land-shuffled run equals the full run to four decimals for all three gases, so the land column is never used by the trees. The small gap for "land off" on CO (12.5 % against 13.1 %) comes from the changed column count altering the random feature sampling, not from land information.

### Land change → CO and SO2 (exploratory; changes no decision)

| | CO | SO2 |
|---|---|---|
| G1 detectable effect, pooled 2019→2023 | **0.13 %/pp, passes** | 4.73 %/pp, fails |
| Own built-up effect, panel [95 % CI] | +0.005 [−0.005, 0.017] | +1.07 [−0.36, 2.51] |
| Long difference | +0.007 [−0.008, 0.022] | +1.14 [−0.77, 2.93] |
| Lead L(t+1), should be 0 | +0.016 [−0.012, 0.040], clean | −0.34 [−2.66, 1.46], uninformative |

- **CO** has a well-powered, clean design: any land effect on CO is below about 0.02 % per pp of built-up gain. The reason is physical: CO columns are spatially smooth and background-dominated, so a 4 km block's land can't move them. This is a precise null, not a failure of power.
- **SO2** can detect only effects above about 4.7 %/pp, so its null says nothing.
- In the CO panel the previous year's land (L_lag) is positive and significant (+0.052 [0.010, 0.085]). It isn't part of the pre-set clean criterion, and it is reported here for completeness.

## 4c. Extensions X1–X5 (rules frozen in `Land Change/EXTENSIONS_PLAN.md` before any result)

All are exploratory or validation-only. None re-opens a frozen decision, and none reads 2025 results or re-scores 2025.

| Step | Result | Verdict |
|---|---|---|
| **X1. Lockdown-free window 2021→2024** | Detectable effect **1.29 %/pp** (signal sd 0.96 pt vs noise 1.7). Own effect −0.056 [−0.13, 0.03]; long difference −0.19 [−0.32, −0.05] raw, −0.11 [−0.26, 0.06] deweathered; lead +0.16 [−0.09, 0.44] raw, but **+0.45 [0.14, 0.71] deweathered** | G1b fails at both 0.5 and 1.0. By the rule, the confirmatory 2021→2025 run is **not executed**. The v1→v2 processor step lies inside this window |
| **X2. Small CNN (validation)** | NO2 +8.4 % [5.2, 11.5] against Ridge +9.4 %; CO +13.6 % [3.6, 24.1] against Ridge +16.3 %. It overfits fast (best epoch 4 and 2) | Not adopted; recorded 2025 scores stand |
| **X3. Interval recalibration** | On 2024, bands fitted on 2023: lead-1 coverage with global quantiles / pixel-scaled quantiles is NO2 0.794 / 0.794, CO 0.859 / 0.859, SO2 0.799 / 0.799 | Not adopted (better for 1 of 3 gases, needed 2). Pixel scaling isn't the issue: CO covered 0.86 in 2024 but 0.74 in 2025, so variance shifts by year |
| **X4. Night-time lights (emission proxy)** | Own effect +0.11 [−0.04, 0.25] % NO2 per 10 % more light; long difference +0.08 [−0.08, 0.23]; lead +0.25 [−0.25, 0.59], clean; **neighbour +0.93 [0.29, 1.52]**. Deweathered is similar. Optimistic detection limit 0.77 % per 10 % light | Own effect inconclusive; a positive neighbour term is a hypothesis for 2026, not a result |
| **X5. ENSO and IOD for long leads** | Mean skill over leads 2–30: NO2 **−3.1 %** [−5.6, −0.3], CO −1.7 % [−5.4, 1.4], SO2 −1.3 % [−2.1, −0.6]. Damped persistence: +0.35 %, +0.14 %, −0.02 % | Worse than climatology; not adopted |

Notes:
- **X4 caveat.** VIIRS annual V2.1 (2019–2021) and V2.2 (2022 on) are joined, and the radiance step at 2021→22 (+0.9 against +0.4 per year before) is a product-version change that year effects absorb. The G1 noise rule assumes land change only grows, which isn't true for lights, so its "signal sd" is not estimable and I report the unadjusted spread as an optimistic bound.
- **Bootstrap CIs against G1 limits.** The panel CIs (about ±0.1–0.2) are tighter than the G1 limits (0.7–1.3) because the panel uses all year pairs while G1 is the conservative long-difference bound with a land-noise deduction. G1 stays the gate.
- **Not run:** the consistent NO2 reprocessing (needs a Copernicus account and ~300 GB).

## 5. Figures

Figures are in `Land Change/results/figures/`:
- `fig_land_layers_g1_lead.png`: G1 detectable effect and the lead and neighbour checks per layer;
- `fig_skill_vs_lead.png`: skill against lead time;
- `fig_fan_2025.png`: Chennai, Coimbatore and Madurai;
- `fig_skill_by_gas_2025.png`: 2025 test skill for NO2, CO and SO2.

`G1_power/fig_g1_*.png` shows the land-change distributions.

## 6. Limitations

- **Coarse resolution.** S5P pixels are ~5 km, and real built-up change is ~1–2 pp per 4 years per block. Dynamic World label noise (~2 pp) is as large as that signal, which attenuates any slope about ×0.3.
- **Changing processor versions.** Every season uses a different TROPOMI version. Site×year effects absorb uniform steps, but a level-dependent step can't be ruled out. The consistent records that exist don't fit this design (checked 13 Oct 2026):
  - The full-mission reprocessing (RPRO v2.4.0, 1 May 2018 – 25 Jul 2022) is Level-2 orbit files only, on the Copernicus Data Space. That means an account and roughly 300 GB for the Jan–Apr seasons 2019–2022 (my estimate: about 700 orbits at about 450 MB each).
  - The S5P-PAL v2.3.1 reprocessing is also Level-2, and it ends in November 2021.
  - The account-free gridded ESA CCI+ product is 0.2° (5× coarser than the blocks) and ends in December 2021.
  - After July 2022 the operational stream is already v2.4.0 or later, consistent with the reprocessing.
- **Lockdown confounding.** The 2020 lockdown sits inside the land window. Excluding 2020 wasn't pre-registered for the panel, so it is reported only as a diagnostic.
- **No CNN.** The deep-learning step wasn't run.
- **Sites.** The 18 sites are Tamil Nadu cities and towns; Bengaluru and Hyderabad lie outside the clipped export.

## 7. Corrections and integrity notes

- **The first 2025 scoring was invalid.** My 2025 Earth Engine download filled pixels outside Tamil Nadu with 0 instead of NaN. That faked 40–60 % drops at coastal and border sites. The method check had missed it because it compared only the −9999 mask.
- **The fix.**
  - The files were repaired into `data/tnaq_extra_v2`.
  - The check now compares the NaN, −9999 and valid classes, and passes.
  - 2025 was re-scored with **unchanged** frozen settings; validation numbers are identical.
  - The invalid record is kept in `Land Change/results/C3_forecast_INVALID_2025_zero_fill/`.
- **Earlier fixes** during the work:
  - a spurious G1 pass on the 1 km proxy (noise mistaken for signal);
  - overlapping Coimbatore and Tiruppur boxes;
  - a floating-point sign issue in the noise estimate.

## 8. Reproduce

See the README section "Phase 3: Land Change". In short:

```bash
uv sync
uv run aqf-land all
uv run aqf-land panel
uv run aqf-land diag
uv run aqf-land --layer dw_prob_trail3 layers
uv run aqf-land --layer dw_prob_trail3 g1
uv run aqf-land --layer dw_prob_trail3 panel
uv run aqf-land bdecision
uv run aqf-land forecast
uv run aqf-land figures
```

The Earth Engine inputs come from `Land Change/scripts/*.py --project airline-portal`.
