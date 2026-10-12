# ADR-001 Phase 1: dual-clock land-use change → NO2 change (data prep + gate G1)

Run on 12 Oct 2026, Windows PC, with the local Tamil Nadu export (`D:\CLG\SEM-7\FINAL PROJECT\Dataset\Raw`, TNAQ).

> **Plan text.** `plans/adr_001_dual_clock_land_emission.md` couldn't be read from this machine, so this phase
> follows the brief in the project thread. These are my assumptions; check them against the ADR:
> - the season is Jan–Apr;
> - the grid is 4×4 blocks of 0.01° (0.04°, 12×12 per city);
> - dL is the built-up fraction change in pp;
> - the assumed effect is **0.5 % NO2 per pp** (placeholder in `configs/land.yaml`);
> - G1 uses the tuning years only.

```bash
uv run aqf-land all        # cache (~40 s first time), layers, provenance, G1; about 1 min
uv run pytest -q "Land Change/tests"
```

## Years

| Use | Years |
|---|---|
| Tuning / G1 | 2019–2023 |
| Old test record (F5, scored once) | 2024. Not used for any fitting or tuning; kept only in the pair table, flagged `held_out_record` |
| Test | **2025**, not downloaded yet |

2018 S2 is TOA rather than SR, so it is never used.

## Data used

- **NO2:** `COPERNICUS/S5P/OFFL/L3_NO2` `tropospheric_NO2_column_number_density`, as 5-day masked means on a 0.01° grid. That's 473 windows (2018-07-10 → 2024-12-30).
- **S2:** `S2_SR_HARMONIZED` (Cloud Score+ ≥ 0.6), as 5-day medians of 6 bands at 0.01° (≈1 km). That's 351 windows for 2019–2024; 2023 has Jan–Apr only and 2024 Jan–Jul.
- **Cities:** nine boxes of 12×12 S5P blocks, 1,114 blocks in total. Chennai is 12×10 because the export ends at 80.36° E, and the lost columns are sea.

## Results

### Processor versions (`results/P1_provenance/`): verified with Earth Engine

The local files don't record the version, so I queried `PROCESSOR_VERSION` from Earth Engine (project
`airline-portal`; `ee_processor_versions.csv`, by month). The Jan–Apr season of each year uses:

| Season | 2019 | 2020 | 2021 | 2022 | 2023 | 2024 | 2025 (test) |
|---|---|---|---|---|---|---|---|
| Processor | 1.2.2 → 1.3.0/1.3.1 | 1.3.2 | 1.4.0 | 2.3.1 | 2.4.0 → 2.5.0 | 2.6.0 | 2.8.0 |

v2.9.1 ("02.09.01") starts in November 2025. December 2025 also has one image labelled `03.00.DV`, which looks like a development build; it's outside the Jan–Apr season.

- **Every season is a different processor, and v1 → v2 (v2.2.0, July 2021) falls between 2021 and 2022.**
- The NO2 data shows no visible step at that pair: the mean across cities is −1.6 %, while other pairs range from −22 % (lockdown) to +16 %. Its level dependence matches every other pair.
- City-year fixed effects absorb any uniform step. G1 is also reported without the 2021→22 pair.

### Gate G1 (`results/G1_power/`): **FAIL** (pooled MDE 0.73 %/pp against the assumed 0.5)

dL comes from **Dynamic World V1** (same-season mode label, 10 m), fetched from Earth Engine for 2019–2025
(`Land Change/data/ee_landfrac/`). NO2 noise comes from the TNAQ export. Tuning years are 2019–2023 only.

| Pooled | cities | pairs | signal sd(dL) | noise sd | sd(dNO2 resid) | N_eff | **MDE %/pp** | needs signal sd ≥ |
|---|---|---|---|---|---|---|---|---|
| 1 km TNAQ proxy, year pairs | 9 | 4,456 | 0 | 32.3 | 7.8 % | 539 | ∞ | 1.9 |
| Dynamic World, year pairs | 9 | 4,456 | 0.68 | 1.24 | 7.8 % | 539 | 1.39 | 1.9 |
| Dynamic World, 2019→2023 | 9 | 1,114 | 2.32 | 1.18 | 6.6 % | 155 | **0.64** | 3.0 |
| Dynamic World, year pairs | 18 | 8,908 | 0 | 1.90 | 7.1 % | 1,328 | ∞ | 1.1 |
| Dynamic World, 2019→2023 | 18 | 2,227 | 1.40 | 1.98 | 6.7 % | 333 | **0.73** | 2.05 |
| *sensitivity:* city-weighted, 2019→2023 | 18 | 2,227 | | | | | *0.45* | |
| *sensitivity:* city-weighted, 2019→2023 | 9 | 1,114 | | | | | *0.58* | |

What this says:
1. **The 1 km proxy is unusable.** Its dL correlates −0.99 with the vegetation change, so it is drought and greenness, not land use. That's why I went to Earth Engine.
2. **Dynamic World shows real, slow urban growth.** Mean built-up change is +0.4 to +0.9 pp per year. Over 2019→2023 the median block gains +1.4 pp and the 95th percentile +7 pp. The fast growers are Tiruppur, Coimbatore, Karur, Hosur, Dindigul and Salem (median +2.5 to +4.6 pp). Built-up gain comes at the expense of vegetation (correlation −0.7 to −0.99), as it should.
3. **Year-to-year pairs can't work.** One year of growth is smaller than the Dynamic World label noise. Only the long difference has signal.
4. **Adding cities doesn't help by itself.** It took me from 9 to 18 cities, and the pooled MDE got worse (0.64 → 0.73). The new boxes mostly don't grow (Thanjavur, Tiruvannamalai, Nagapattinam: signal 0), and their label noise raises the pooled noise estimate. Weighting cities by their own precision gives 0.45 (18 cities), which would pass. That is a sensitivity only: the gate as defined is the pooled one, and switching estimators after seeing the result would be a forking path. **Decide in the ADR which estimator is G1, before anything about the slope is looked at.**
5. **Chennai alone is hopeless.** Its MDE is 3.0 %/pp: strong spatial correlation (ρ 0.82) leaves N_eff = 10.
6. The assumed effect of 0.5 %/pp is a placeholder. With 1.0 %/pp, the 9- and 18-city long differences both pass.

Box fix: Coimbatore's and Tiruppur's boxes overlapped in the first run (their centres are only 0.37° apart). Tiruppur now
sits 12 km east, and a test enforces that no two boxes overlap. The 9-city rows above are from the first run.

## Decisions of 12 Oct 2026 and what followed

Naresh decided:
- **G1** is the pooled long difference 2019→2023, gated at 0.5 %/pp and also reported at 1.0 %/pp;
- sites are added only by Dynamic World growth;
- if G1 fails, the C2 panel runs anyway, labelled exploratory.

### Site selection (`SITE_SELECTION.md`, committed before any slope)

- **Result: 0 of 10 candidate tiles qualify.** Their growth scores run from −0.25 to +0.84 pp, against the frozen 1.5 pp threshold.
- The rest of Tamil Nadu, outside the 18 boxes, isn't growing fast enough to add power, so no sites were added.
- I didn't loosen the threshold after seeing this (`results/G1_power/site_selection.csv`).

### G1, final

| Design | MDE (80 % power) | at 0.5 %/pp | at 1.0 %/pp |
|---|---|---|---|
| Pooled long difference 2019→2023, 18 sites | **0.73 %/pp** | **FAIL** | PASS |

### ERA5 (done, through Earth Engine; no CDS account needed)

- **Source:** `scripts/ee_era5.py` reads `ECMWF/ERA5/HOURLY`, the same ERA5 single levels as CDS.
- **Reduction:** 06–09 UTC hours, each site box widened by 0.25°. One row per site and 5-day window, 2019–2025.
- **Output:** `results/P1_era5/era5_city_windows.csv`, with the same columns as `cds_era5.py`.
- **Not used:** 2024 and 2025 aren't used by any model.

### C2 panel, EXPLORATORY (`results/C2_panel_EXPLORATORY/`)

It uses tuning years 2019–2023 only. The model has pixel and site×year fixed effects plus own and neighbour built-up (SLX). Uncertainty comes from a spatial-block bootstrap (4×4 pixels, 500 draws), alongside a spatial placebo, a lead/lag check and the long difference. I ran it on raw NO2 and on ERA5-deweathered NO2.

| Outcome | own built-up b (%/pp) [95 % CI] | long difference b | lead L(t+1), must be ≈ 0 | neighbour WL |
|---|---|---|---|---|
| raw | −0.004 [−0.08, 0.09] | −0.03 [−0.15, 0.10] | **+0.31 [0.14, 0.45]** | +0.27 [0.05, 0.48] |
| deweathered | −0.002 [−0.08, 0.07] | +0.00 [−0.13, 0.14] | **+0.35 [0.19, 0.46]** | +0.41 [0.20, 0.60] |

Reading:
1. **No own-pixel effect is visible.** The bootstrap CIs are narrower than the G1 MDE (0.73), because the MDE is conservative: it counts every negative dL as noise.
2. **Attenuation.** Dynamic World label noise (≈ 2 pp) is larger than the real 4-year change (≈ 1.4 pp), so the estimates are attenuated toward zero by roughly ×0.3. An attenuation-corrected upper bound is about 0.3 %/pp; this is a rough figure.
3. **The lead test fails in both versions.** Next year's built-up "predicts" this year's NO2, and the neighbour term is significant. Neither can be causal, so **the panel is confounded and none of its land coefficients should be read as effects.**
   - A likely cause is that Dynamic World labels respond to the same conditions as NO2, for example haze or dry-season surface brightness.
   - Another is that growth areas differ in unmodelled ways.
4. Deweathering changes nothing important. With site×year fixed effects, site-wide weather is already absorbed.

## End-to-end run, 12 Oct 2026 evening

**Summary.** The consolidated write-up is [`docs/FINAL_REPORT.md`](../docs/FINAL_REPORT.md), with a chapter draft in
[`docs/THESIS_CHAPTER_DRAFT.md`](../docs/THESIS_CHAPTER_DRAFT.md). Each step's rules were committed before its results.

| Step | Rule file | Outcome |
|---|---|---|
| Lead-test diagnostics | `DIAGNOSTICS_PLAN.md` | Lead driven by the **2020 lockdown** outcome year (0.31 → 0.01 without it); aerosol (AOD) flagged by D1 (2019–2021 covariates; 2022–2024 timed out in EE); not edges, single sites or shared spatial residuals (`results/D_lead_diagnostics/`) |
| Steadier layers | `DIAGNOSTICS_PLAN.md` | G1 MDE: DW label 0.73, DW prob. trailing-3y 1.52, Open Buildings 4.17 %/pp. Only DW prob. is lead-clean on raw NO2; with deweathered NO2 its lead fails too (`results/layers/`) |
| Land term B | `DECISION_RULES.md` §B | **B = 0**: G1 fails on every layer (`results/B_decision.json`) |
| Forecast (validation) | `DECISION_RULES.md` §C–D | Lead 1: Ridge **+9.4 %** [6.4, 13.4], LightGBM +6.5 %, damped persistence +5.4 %. ERA5 helps slightly; land adds nothing. CNN not run (PyTorch download stopped) |
| 2025 test (frozen) | `DECISION_RULES.md` §E | Lead 1 **+3.6 %** [−1.3, 6.8]; leads 2–30 ≈ 0; 10–90 % coverage 0.81 (`results/C3_forecast/`) |
| Integrity | – | The first 2025 scoring is invalid: the EE download filled outside-AOI pixels with 0. Files were repaired into `data/tnaq_extra_v2`, the check was fixed, and 2025 was re-scored unchanged. The record is kept in `results/C3_forecast_INVALID_2025_zero_fill/` |

## Still to do (needs your accounts)

| Step | Script | Needs |
|---|---|---|
| ~~Exact processor versions~~ | **done** (`ee_s5p_versions.py --project airline-portal`) | |
| ~~Real dL~~ | **done**: `ee_s2_landfrac.py --project airline-portal --years 2019 … 2025 --classifier dynamicworld --direct` | |
| ~~ERA5~~ | **done via Earth Engine** (`ee_era5.py`). `cds_era5.py` stays as the CDS alternative. | |
| ~~2025 test-year NO2~~ | **done**: `ee_s5p_no2_windows.py --project airline-portal --year 2025`. 73 windows went to `data/tnaq_extra/` (22 MB), not `Dataset/`. Same method as the TNAQ export: `--check` on 2021-07-04 and 2024-01-05 reproduced the local files bit for bit (same mask, max diff 0). 2025 is loaded but not used by any step yet; it's for the final scoring only. | |

The ERA5 reducer was tested end-to-end on a synthetic NetCDF only; no ERA5 data exists on this machine.

## Files

- `src/aqf_land/tnaq.py`: city boxes, cached sub-pixel reads, block means.
- `src/aqf_land/layers.py`: seasonal composite, classification, class fractions, seasonal NO2, the pixel-year pair table, and the EE CSV loader.
- `src/aqf_land/provenance.py`, `g1.py`, `era5.py`, `cli.py`.
- `scripts/`: Earth Engine and CDS jobs.
- `results/P1_layers/` is gitignored (rebuildable); `P1_provenance/` and `G1_power/` are tracked.
