# aqf: satellite air-pollution analysis and forecasting over Chennai

Sentinel-5P NO2/CO/SO2 column maps (5-day composites, 2019–2024, ~5 km) are the pollution signal.
Sentinel-2 surface reflectance (100 m) is land-surface context. Raw data lives in `Dataset/` (see `Dataset/dataset.md`).

Current status, key results and next steps: [`PROJECT_STATUS.md`](PROJECT_STATUS.md).

The project is a uv workspace. The shared core package `aqf` (`src/aqf/`) holds the config, data loading and
data-quality rules. Each phase lives in its own folder with its own code, tests, results and report.

## Setup

```bash
uv sync              # Python 3.13; installs the core and every phase package, pinned in uv.lock
uv run pytest -q     # tests of the core and all phases
```

## Phase 1: Data Analysis (`Data Analysis/`)

```bash
uv run aqf-analysis all          # all steps below, about 1 minute; add --no-hash to skip SHA-256
uv run aqf-analysis inventory    # 01 manifest/file integrity, grids, provenance, S2 availability
uv run aqf-analysis cube         # build data/processed/s5p_cube.npz (438 x 3 x 13 x 13)
uv run aqf-analysis quality      # 02 coverage, zeros, unusable windows
uv run aqf-analysis trends       # 03 time series, seasonal cycle, STL, Mann-Kendall/Sen trends, lockdown
uv run aqf-analysis spatial      # 04 mean/seasonal maps, hotspots, industrial hubs
uv run aqf-analysis s2           # 05 Sentinel-2 coverage and spectral indices (streamed)
```

The findings are in [`Data Analysis/REPORT.md`](<Data Analysis/REPORT.md>). Figures and tables go to
`Data Analysis/results/<NN_step>/`.

## Phase 2: Forecasting (`Forecasting/`)

Task A forecasts the next 5-day map. Task B is a 150-day (30-window) outlook. Both are for NO2, CO and SO2 on the 12×12 grid.
The plan is [`Forecasting/PLAN.md`](Forecasting/PLAN.md) and the findings so far are in [`Forecasting/REPORT.md`](Forecasting/REPORT.md).

```bash
uv run aqf-forecast all              # F0 + F1 + F2 + F4, about 4 min
uv run aqf-forecast forecastability  # F0 autocorrelation, cross-gas, neighbours, S2 vs long-term mean
uv run aqf-forecast baselines        # F1 climatology / persistence / damped persistence, tasks A and B, validation 2023
uv run aqf-forecast learned          # F2 LightGBM + Ridge vs baselines, about 3 min
uv run aqf-forecast ablations        # F4 ablations, rolling-origin 2021-2023, NE-monsoon gate (after F2)
```

## Layout

```
configs/forecast.yaml                  horizon, lookback, quantiles, bootstrap, forecastability settings
configs/data.yaml                      paths, gases, display units, seasons, splits, quality rules, hub coordinates
src/aqf/data.py                        core: config, manifest, raster reading, S5P cube
src/aqf/series.py                      core: zero handling, per-window quality table, area-mean and monthly series
src/aqf/anomaly.py                     core: per-pixel harmonic climatology (fit on train only)
src/aqf/s2_static.py                   core: static S2 features aggregated to the S5P grid (train years)
tests/                                 core tests
Data Analysis/
  REPORT.md                            written findings of phase 1
  src/aqf_analysis/                    cli, stats (Mann-Kendall, Sen, FDR), plotting, one module per step
  tests/                               analysis tests
  results/                             generated figures and tables (gitignored)
Forecasting/
  PLAN.md, REPORT.md                   phases F0-F5 and findings
  src/aqf_forecast/                    problem, baselines, features, learned (F2), ablations (F4), metrics, forecastability, evaluate
  tests/                               climatology, no-lookahead, metric and alignment tests
  results/                             generated figures and tables (gitignored)
```
