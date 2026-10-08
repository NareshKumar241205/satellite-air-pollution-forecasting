# aqf: satellite air-pollution analysis and forecasting over Chennai

Sentinel-5P NO2/CO/SO2 column maps (5-day composites, 2019–2024, ~5 km) are the pollution signal.
Sentinel-2 surface reflectance (100 m) is land-surface context. Raw data lives in `Dataset/` (see `Dataset/dataset.md`).

## Setup

```bash
uv sync              # Python 3.13, dependencies pinned in uv.lock
uv run pytest -q
```

## Data analysis

```bash
uv run aqf analysis          # all steps below, about 1 minute; add --no-hash to skip SHA-256
uv run aqf inventory         # 01 manifest/file integrity, grids, provenance, S2 availability
uv run aqf cube              # build data/processed/s5p_cube.npz (438 x 3 x 13 x 13)
uv run aqf quality           # 02 coverage, zeros, unusable windows
uv run aqf trends            # 03 time series, seasonal cycle, STL, Mann-Kendall/Sen trends, lockdown
uv run aqf spatial           # 04 mean/seasonal maps, hotspots, industrial hubs
uv run aqf s2                # 05 Sentinel-2 coverage and spectral indices (streamed)
```

Everything goes to the `Data Analysis/` folder, with one sub-folder per step (`01_inventory/` … `05_sentinel2/`)
holding CSV tables and PNG figures. The findings are written up in
[`Data Analysis/REPORT.md`](<Data Analysis/REPORT.md>).

## Layout

```
configs/data.yaml        paths, gases, display units, seasons, splits, quality rules, hub coordinates
src/aqf/data.py          config, manifest, raster reading, S5P cube
src/aqf/series.py        zero handling, per-window quality table, area-mean and monthly series
src/aqf/stats.py         seasonal Mann-Kendall, seasonal Sen slope, Benjamini-Hochberg FDR
src/aqf/plotting.py      shared figure style, gas maps with hubs
src/aqf/analysis/        inventory, quality, trends, spatial, s2 (one module per step)
Data Analysis/           analysis results (generated) + REPORT.md (written findings)
```
