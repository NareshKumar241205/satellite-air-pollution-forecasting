# aqf: satellite air-pollution analysis and forecasting over Chennai

Sentinel-5P NO2/CO/SO2 column maps (5-day composites, 2019–2024, ~5 km) are the pollution signal.
Sentinel-2 surface reflectance (100 m) is land-surface context. Raw data lives in `Dataset/` (see `Dataset/dataset.md`).

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

## Layout

```
configs/data.yaml                      paths, gases, display units, seasons, splits, quality rules, hub coordinates
src/aqf/data.py                        core: config, manifest, raster reading, S5P cube
src/aqf/series.py                      core: zero handling, per-window quality table, area-mean and monthly series
tests/                                 core tests
Data Analysis/
  REPORT.md                            written findings of phase 1
  src/aqf_analysis/                    cli, stats (Mann-Kendall, Sen, FDR), plotting, one module per step
  tests/                               analysis tests
  results/                             generated figures and tables (gitignored)
```
