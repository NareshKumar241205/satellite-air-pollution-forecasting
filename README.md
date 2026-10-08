# aqf: satellite air-pollution analysis and forecasting over Chennai

Sentinel-5P NO2/CO/SO2 column maps (5-day composites, 2019–2024, ~5 km) are the pollution signal.
Sentinel-2 surface reflectance (100 m) is land-surface context. Raw data lives in `Dataset/` (see `Dataset/dataset.md`).

Current status, key results and next steps: [`PROJECT_STATUS.md`](PROJECT_STATUS.md).

## Results at a glance

**Data analysis** ([`Data Analysis/REPORT.md`](<Data Analysis/REPORT.md>)):
- Exact zeros in the S5P maps are missing retrievals, not clean air.
- No gas has a significant 2019–2024 trend.
- NO2 fell 35.5 % during the 2020 lockdown.
- The NO2 hotspot is on the Manali–Ambattur axis, and SO2 peaks toward Ennore.

**Forecasting** ([`Forecasting/REPORT.md`](Forecasting/REPORT.md)): 5-day and 150-day forecasts of the 12×12 maps.
The models were chosen with rolling-origin validation over 2021–2023 and scored **once** on the held-out year 2024:

| Gas | 5 days | 10–30 days | 35–150 days | 80 % interval coverage (2024) |
|---|---|---|---|---|
| NO2 | damped persistence +3.1 % [−0.4, 6.1] | Ridge −0.5 % | seasonal climatology | 0.82 |
| CO | damped persistence +3.8 % [−7.6, 10.1] | climatology | climatology | 0.85–0.87 |
| SO2 | climatology | climatology | climatology | 0.81 |

The numbers are MSE skill against the seasonal climatology, with 95 % CIs in brackets. Beyond the seasonal cycle, predictability is small and
year-dependent, so the dependable long-range forecast is the seasonal climatology with calibrated bands. The **Jan–May 2025
outlook** (GeoTIFF per gas, CSVs, figures) is in [`Forecasting/results/F5_final/outlook_2025/`](Forecasting/results/F5_final/outlook_2025/).

## Data (not included in this repository)

The raw dataset (1.9 GB) is not in git. Place it at `Dataset/data/` with this layout (described in
[`docs/dataset.md`](docs/dataset.md)):

```
Dataset/data/
  dataset_manifest.csv          438 five-day windows, 2019-01-01 .. 2024-12-25
  s5p_composites/s5p_YYYYMMDD.tif   Sentinel-5P, 3 bands (NO2, CO, SO2; mol/m^2), 13x13 at ~5 km, EPSG:4326
  s2_composites/s2_YYYYMMDD.tif     Sentinel-2, 12 bands (B2-B12, NDVI, NDBI, NDMI), 558x558 at 100 m (290 files)
```

Without the data, `uv run pytest` still runs the unit tests and skips the tests that need it.

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

Task A forecasts the next 5-day map. Task B is a 150-day (30-window) outlook. The **Jan–May 2025 outlook** (GeoTIFF/CSV/figures)
is in `Forecasting/results/F5_final/outlook_2025/`. Both are for NO2, CO and SO2 on the 12×12 grid.
The plan is [`Forecasting/PLAN.md`](Forecasting/PLAN.md) and the findings so far are in [`Forecasting/REPORT.md`](Forecasting/REPORT.md).

```bash
uv run aqf-forecast all              # F0 + F1 + F2 + F4, about 4 min
uv run aqf-forecast forecastability  # F0 autocorrelation, cross-gas, neighbours, S2 vs long-term mean
uv run aqf-forecast baselines        # F1 climatology / persistence / damped persistence, tasks A and B, validation 2023
uv run aqf-forecast learned          # F2 LightGBM + Ridge vs baselines, about 3 min
uv run aqf-forecast ablations        # F4 ablations, rolling-origin 2021-2023, NE-monsoon gate (after F2)
uv run aqf-forecast final            # F5: score 2024 ONCE (locked afterwards) + Jan-May 2025 outlook
uv run aqf-forecast outlook          # F5 outlook only
```

## Viewing any GeoTIFF

Open [`tools/view_tiff.ipynb`](tools/view_tiff.ipynb) in VS Code or Jupyter and choose the kernel **Python (Scratch .venv)**
(*Select Kernel → Jupyter Kernel… → Python (Scratch .venv)*). That kernel is registered once per machine with:

```bash
uv run python -m ipykernel install --user --name aqf-scratch --display-name "Python (Scratch .venv)"
```

If the kernel list is empty in Code - OSS, see [`tools/README.md`](tools/README.md) (proposed-API fix + a JupyterLab fallback:
`uv run --with jupyterlab jupyter lab tools/view_tiff.ipynb`).
Set `IMAGE_PATH` in the settings cell (relative to `Scratch/`), then *Run All*. The other settings are `BANDS`, `RGB`
(Sentinel-2 true colour), `VMIN`/`VMAX` and `SAVE_TO`. The last cell lists every `.tif` in the project.

```bash
uv run jupyter lab tools/view_tiff.ipynb      # only if you don't use VS Code; needs `uv add --dev jupyterlab`
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
tools/view_tiff.ipynb                  notebook to view any GeoTIFF (set IMAGE_PATH)
Forecasting/
  PLAN.md, REPORT.md                   phases F0-F5 and findings
  src/aqf_forecast/                    problem, baselines, features, learned (F2), ablations (F4), final (F5 test + outlook), metrics, forecastability, evaluate
  tests/                               climatology, no-lookahead, metric and alignment tests
  results/                             generated figures and tables (gitignored)
```
