# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A final-year project on satellite air-pollution forecasting over Chennai. It was rebuilt from scratch on 2026-10-08,
and the earlier code isn't part of this repo. Phase 1 (data analysis) is done (`Data Analysis/REPORT.md`). Phase 2 (forecasting) follows `Forecasting/PLAN.md`:
F0 and F1 are done (`Forecasting/REPORT.md`) and F2 (LightGBM/Ridge) is next. Read both reports before modelling work.

## Commands

```bash
uv sync                                        # installs the core `aqf` + workspace members (Data Analysis)
uv run pytest -q                               # 14 tests, about 5 s; tests/test_data.py skips if Dataset/ is absent
uv run pytest "Data Analysis/tests/test_stats.py::test_fdr_bh_matches_hand_computation" -q
uv run aqf-analysis all --no-hash              # every analysis step, about 1 min
uv run aqf-analysis <inventory|cube|quality|trends|spatial|s2>
uv run aqf-forecast <forecastability|baselines|all>   # Forecasting F0/F1, about 40 s
```

## Architecture

- **uv workspace.** The root package `aqf` (`src/aqf/`) is the shared core that every phase reuses. Each phase is a
  workspace member in its own folder holding its code, tests, results and report. `Data Analysis/` is the package
  `aqf-analysis` (module `aqf_analysis`). A new phase (forecasting, for example) becomes a new member folder that is added to
  `[tool.uv.workspace] members` and to the root `dev` dependency group, and imports from `aqf`.
- `configs/data.yaml` holds all paths (relative to this directory) and every data choice. Change the config rather than
  hard-coding values.
- `aqf/data.py::build_s5p_cube` stacks all 438 S5P rasters into `data/processed/s5p_cube.npz` (NaN = nodata) after asserting
  that every file shares one grid. Sentinel-2 is never stacked: at 558×558×12 per window it doesn't fit in this machine's
  ~7 GB of RAM, so `aqf_analysis/s2.py` streams one file at a time.
- `aqf/series.py` is the single place where data-quality rules apply. `clean_values` turns zeros into NaN,
  `window_table` flags usable windows (≥ `min_coverage` valid pixels), and `analysis_values` gives the cleaned and masked
  cube that every pixel-level analysis must use. Don't read `cube.values` directly in an analysis.
- **Forecasting** (`Forecasting/src/aqf_forecast/`): `problem.py` loads the cleaned cube and splits. `fit()` fits the
  climatology (`aqf/anomaly.py`) and the damped-persistence coefficients on a window mask. Every model returns
  predicted values of shape (T origins, H leads, G, 13, 13) for all origins at once. A pair's split is that of its
  *target* window (`target_obs`, `_target_split`). Intervals and CRPS come from leave-one-train-year-out residuals.
  `tests/test_forecast.py::test_no_lookahead_in_baselines` guards against future leakage, so extend it for every new model.
- Each `aqf_analysis/<step>.py` has a `run()` that writes to `Data Analysis/results/<NN_step>/` (`paths.outputs`).
  `aqf_analysis/cli.py` wires them up.

## Data facts that bite

- S5P: the 13×13 grid has 12×12 valid pixels. Row 12 and column 0 are always `-inf`.
- **Exact zeros are missing retrievals, not clean air** (evidence in `Data Analysis/REPORT.md` §2). SO2 is 39 % zeros,
  and only 79 % of its windows are usable.
- Sentinel-2 exists for only 290 of 438 windows (monsoon cloud gaps). About 40 % of each raster is sea. NDMI is exactly −NDBI.
- GDAL logs a harmless photometric/ExtraSamples warning for the S2 TIFFs. The CLI silences the `rasterio` logger.

## Rules

- 2024 (test) is scored once, in F5, with frozen configs. Everything before that is evaluated on 2023 only.
- `Dataset/` is raw and immutable. `data/`, `Data Analysis/results/` and `Forecasting/results/` are derived and gitignored. Never commit them.
- Split by year: train 2019–2022, val 2023, test 2024. Anything fitted (scalers, climatology, S2 composites) uses train
  years only. Descriptive analysis may use all years.
