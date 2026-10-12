# ADR-001 plan and status (updated 13 Oct 2026)

Branch `adr-001-phase1`. The original plan text (`plans/adr_001_dual_clock_land_emission.md`) is in the project's shared files and was never readable from this PC, so this file records the plan as it was actually executed from the coordinator's brief and Naresh's decisions. It supersedes earlier status notes in `REPORT.md`.

## 1. Status at a glance

| Phase | Item | Status | Where |
|---|---|---|---|
| A | Inventory and cleanup | Done. Nothing needed trashing (fresh clone) | `CLEANUP.md` |
| P1 | Data prep: TNAQ city boxes, yearly land fractions, same-season NO2 | Done, 18 sites | `src/aqf_land/` |
| P1 | TROPOMI processor versions | Done (Earth Engine, 2019–2025) | `results/P1_provenance/` |
| P1 | ERA5 scalars per site and window | Done, via Earth Engine (no CDS account) | `results/P1_era5/` |
| G1 | Power gate | **Failed** at 0.5 %/pp (0.73); passes at 1.0 | `results/G1_power/` |
| Sites | Growth-only site rule | Done: 0 of 10 candidates qualify | `SITE_SELECTION.md` |
| C2 | Exploratory panel + lead test | Done, confounded by the 2020 lockdown | `results/C2_panel_EXPLORATORY/` |
| D1–D7 | Lead-test diagnostics | Done (D7 post hoc) | `DIAGNOSTICS_PLAN.md` |
| Layers | Steadier land layers | Done (DW probability, Open Buildings) | `results/layers/` |
| B | Land term in the forecast | **B = 0** (rule applied mechanically) | `results/B_decision.json` |
| C3–E | NO2 forecast ladder, fusion, 2025 test | Done, 2025 scored once | `results/C3_forecast/` |
| G | CO and SO2 (land check, forecast, 2025 test) | Done, each scored once | `results/gas/` |
| F | Reports, chapter draft, README, deck, figures | Done | `docs/`, `results/figures/` |
| X1–X5 | Lockdown-free window, CNN check, band recalibration, night lights, climate indices | Done; none adopted (see `FINAL_REPORT.md` §4c) | `EXTENSIONS_PLAN.md`, `results/` |
| Push | Branch pushed, pull request into `main` | Branch pushed; **pull request not opened** (no `gh` login) | see §4 |

## 2. Frozen decisions (don't change after seeing results)

- **G1 estimator:** pooled long difference 2019→2023. The city-weighted version is a sensitivity only.
- **Effect to power for:** 0.5 %/pp, also reported at 1.0.
- **Land term:** B = 0 for NO2, CO and SO2. Re-opening it would need new data, not new estimators.
- **Sites:** the 18 in `configs/land.yaml`. The site rule found nothing to add.
- **Years:**
  - land analysis 2019–2023;
  - forecast train 2019–2022, validation 2023–2024;
  - **2025 is spent**: it has been scored once for each of NO2, CO and SO2. Any further modelling must not use 2025 for tuning, and any new test needs new data (2026).
- **Model choice (lead 1):** the best validation MSE in the pre-registered ladder. Ridge's α = 1000 is at the edge of the grid; the grid wasn't widened.
- **Records kept, not hidden:** the invalid first 2025 scoring (`C3_forecast_INVALID_2025_zero_fill/`), the spurious early G1 pass, and the withdrawn aerosol claim.

## 3. Results to quote

| | NO2 | CO | SO2 |
|---|---|---|---|
| Land detection limit (G1) | 0.73 %/pp | 0.13 %/pp | 4.73 %/pp |
| Land effect | none above the limit | below ≈ 0.02 %/pp | untestable |
| Validation skill, lead 1 | +9.4 % [6.4, 13.4] | +16.3 % [5.2, 27.0] | ≈ 0 |
| **2025 test, lead 1** | **+3.6 % [−1.3, 6.8]** | **+25.8 % [12.0, 38.6]** | +0.1 % [0.0, 0.2] |
| Leads ≥ 10 days | ≈ 0 | ≈ 0 | ≈ 0 |
| 10–90 % coverage, 2025 | 0.81 | 0.74 | 0.94 |

## 4. What is left, in order

**Needs Naresh (nothing else can proceed without these):**
1. **Open the pull request** into `main`, then review the diff. Link: https://github.com/NareshKumar241205/satellite-air-pollution-forecasting/pull/new/adr-001-phase1
2. **Fill the deck cover placeholders** (name, register number, guide, date).
3. **Read the chapter draft** (`docs/THESIS_CHAPTER_DRAFT.md`) and tell me what your guide wants changed in framing. The honest headline is "no land→NO2 effect above 0.73 %/pp is detectable; NO2 is forecastable only in the next window, CO more so".

**Optional improvements.** Done (12 Oct–13 Oct): the lockdown-free window, the CNN check, band recalibration, night lights and climate indices, none adopted. Still open:

**Original optional list, for reference (each a new, pre-registered step, never tuned on 2025):**

| Item | Effort | What it would change |
|---|---|---|
| Land window avoiding the lockdown (2021→2024, using 2025 only as the final check) | 1–2 h | Removes the 2020 confounder by design instead of post hoc; may change G1 |
| Recalibrate the interval bands (CO coverage 0.74, SO2 0.94) using validation quantiles per gas and lead | ~1 h | Bands only; no change to point forecasts |
| Install PyTorch and run the small CNN | 1–2 h plus the download | May beat Ridge on validation; the rule allows it |
| Emission-weighted land variables (roads, industrial land use) | days | The only route likely to give detectable land effects at 5 km |
| Consistently reprocessed NO2 (v2.4.0 Level-2) | days, ~300 GB, Copernicus account | Removes the processor-version caveat |
| Climate indices (IOD, ENSO) as long-lead predictors | ~half a day | Possible skill beyond 10 days; currently zero |

**Housekeeping:** ignored folders `data/` and `Land Change/results/layers/*` are rebuildable; the 2025 downloads are in `data/tnaq_extra_v2` (not in git). `Dataset/` was never modified.

## 5. Reproduce

See the README section "Phase 3: Land Change". Each gas uses `--gas CO` or `--gas SO2`; each layer uses `--layer dw_prob_trail3` or `--layer obt`.
