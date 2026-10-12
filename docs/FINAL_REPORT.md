# Final report: land-use change and Sentinel-5P NO2 over Tamil Nadu, with a 5-day to 150-day forecast

Branch `adr-001-phase1`, 12 Oct 2026. Every number below comes from files in `Land Change/results/`. Decision rules were committed before the results they govern:
- `Land Change/SITE_SELECTION.md`
- `Land Change/DIAGNOSTICS_PLAN.md`
- `Land Change/DECISION_RULES.md`

## 1. Question and answer

**Question.** Does land-use change (built-up growth) measurably change tropospheric NO2 in Sentinel-5P? Can a
satellite NO2 forecast for Tamil Nadu cities use land change, and how far ahead is it skilful?

**Answer.**
1. **No land→NO2 effect is detectable.** On the best available design, the smallest effect that could have been detected (80 % power) is **0.73 % NO2 per percentage point of built-up gain**. Every estimate is near zero; for example, own built-up gives −0.004 %/pp [−0.08, 0.09]. The design checks fail, though: the lead placebo is non-zero, driven by the 2020 lockdown year and correlated with aerosol. So the null can't be read as "no effect" either. The thesis result is the detection limit.
2. **The land term in the forecast is therefore zero (B = 0),** under a rule fixed in advance.
3. **Forecasts beat climatology only for the next 5-day window.**
   - Validation 2023–24: **+9.4 %** MSE skill [6.4, 13.4] (Ridge).
   - Test 2025, scored with frozen settings: **+3.6 %** [−1.3, 6.8].
   - From 10 days to 150 days, skill is ≈ 0, so the dependable long-range forecast is the per-pixel seasonal climatology with calibrated bands. The 10–90 % band covered 81 % of 2025 pixel values.

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
- **The lockdown year.** Dropping outcome year 2020 removes the lead: 0.31 → **0.01** (−97 %). Jan–Apr 2020 contains the COVID lockdown, whose NO2 drop varied across space. Pixels that later gained built-up area had different lockdown drops.
- **Aerosol (D1, covariates for 2019–2021 only).** Within site and year, AOD correlates 0.61 with the DW built-up label, and next-year AOD correlates with both L(t+1) (0.63) and the NO2 residual (−0.12). That crosses the pre-set |r| ≥ 0.1 flag. Urban, hazy pixels are also built-up pixels, so a common aerosol factor can't be excluded.

**Ruled out:**
- **Single sites:** leave-one-site-out moves the lead by at most ±24 %.
- **Box edges:** dropping the outer ring changes it by −11 %.
- **Shared spatial structure:** the neighbour cross-correlation of the residuals is 0.02.

**Partly explained:**
- **Pixel trends:** pixel-specific trends shrink the lead to 0.11.
- **Level dependence:** yearly label change is larger in high-NO2 pixels (+0.8 to +2.0 pp per 100 log-units in 2019–2022), which is consistent with growth on the urban fringe.

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
- **Skill by lead:** damped persistence falls to +1.7 % at 10 days and below +1 % from 15 days.
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
- **Bias:** 2025 NO2 ran a median 7 % below climatology.

## 5. Figures

Figures are in `Land Change/results/figures/`:
- `fig_land_layers_g1_lead.png`: G1 detectable effect and the lead and neighbour checks per layer;
- `fig_skill_vs_lead.png`: skill against lead time;
- `fig_fan_2025.png`: Chennai, Coimbatore and Madurai.

`G1_power/fig_g1_*.png` shows the land-change distributions.

## 6. Limitations

- **Coarse resolution.** S5P pixels are ~5 km, and real built-up change is ~1–2 pp per 4 years per block. Dynamic World label noise (~2 pp) is as large as that signal, which attenuates any slope about ×0.3.
- **Changing processor versions.** Every season uses a different TROPOMI version. Site×year effects absorb uniform steps, but a level-dependent step can't be ruled out. A consistent reprocessing (S5P-PAL) wasn't available in Earth Engine.
- **Lockdown confounding.** The 2020 lockdown sits inside the land window. Excluding 2020 wasn't pre-registered for the panel, so it is reported only as a diagnostic.
- **Incomplete aerosol check.** D1 covers 2019–2021 only: Earth Engine timed out on 2022–2024.
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
