# 01 — Data analysis: Chennai Sentinel-5P / Sentinel-2, 2019–2024

To reproduce, run `uv run aqf analysis` (about 1 minute). It writes every number and figure below to
`outputs/analysis/<step>/`. The analysis is descriptive and covers all six years. Nothing in it is fitted
for a model, except the Sentinel-2 long-term maps, which use only the train years 2019–2022.

## 1. What the data actually is (step 01_inventory)

| | Sentinel-5P (pollution) | Sentinel-2 (land surface) |
|---|---|---|
| Files | 438, one per 5-day window, none missing | **290** (`dataset.md` says 438; 148 windows are `s2_cloud_gap`) |
| Grid | 13×13 at 0.0449° (~5 km), EPSG:4326 | 558×558 at 0.0009° (~100 m) |
| Area of interest | **12×12 = 144 pixels**: row 12 and column 0 are always nodata (`-inf`) | ~40 % of the raster is sea (Bay of Bengal) |
| Bands | NO2, CO, SO2 column density, mol/m² | B2–B8, B11, B12, NDVI, NDBI, NDMI |

- The windows form a continuous 5-day grid from 2019-01-01 to 2024-12-25 (73 per year; 2020 has 74 and 2024 has 72). There are no duplicates.
- All S5P files share one grid, and so do all S2 files. SHA-256 hashes of every file are in `01_inventory/provenance.csv`.
- **Sentinel-2 is often missing in the monsoon.** Only 31 % of July windows and 44 % of August windows have an image, against 95–97 % in Feb–Apr (`fig01_availability.png`). Any model that needs S2 at every time step will have gaps during the SW monsoon.
- **NDMI is exactly −NDBI.** Both are computed from B8 and B11 with opposite signs, so one of the two bands is redundant.
- `dataset.md` also describes S5P as "~11×11". It is 13×13.

## 2. Data quality: zeros are missing data (step 02_quality)

The S5P maps contain many **exact zeros**: 7 % of NO2 pixels, 11 % of CO pixels and 39 % of SO2 pixels. The zeros are missing
retrievals (cloud-screened pixels written as 0), not clean air:
- **CO:** the smallest non-zero value is 0.0145 mol/m² and the median is 0.038. No values lie between 0 and 0.0145, and a
  total CO column is never physically zero.
- 76–79 % of NO2 and CO zeros sit in spatially clustered patches.
- Zeros peak in the cloudy monsoon months (CO: 23–24 % in Jul–Aug, 2–3 % in Feb–Apr) (`fig03`).
- 2020 has the most zeros of any year for all three gases (CO 18 %, NO2 12 %, SO2 45 %).

**Decision:** zeros are treated as missing (`quality.zeros_as_missing` in `configs/data.yaml`). A gas map is
usable when at least 50 % of the 144 pixels are valid. That leaves this many usable windows:

| | NO2 | CO | SO2 |
|---|---|---|---|
| Usable windows (of 438) | 411 (94 %) | 409 (93 %) | **344 (79 %)** |
| Entirely-zero maps | 5 | 11 | 9 |

Counting zeros as values would bias every mean down, most of all in the monsoon and in 2020. In the first
pass it inflated the SO2 lockdown drop from −21 % to −38 %. SO2 also contains retrieval noise that
was floored at zero, so its means are somewhat biased upward. **SO2 is the least reliable of the three gases.**

## 3. Temporal behaviour and trends (step 03_trends)

### Seasonality dominates; long-term trends are not detectable

Variance shares come from STL decomposition of the 5-day series with a one-year period (`stl_strength.csv`, `fig10`):

| Gas | Seasonal share | Trend share | Residual share | Peak month | Low month | Seasonal amplitude |
|---|---|---|---|---|---|---|
| NO2 | 0.55 | 0.05 | 0.62 | Jun | Apr | 36 % of the mean |
| CO | **0.83** | 0.01 | 0.20 | Mar | Jul | 44 % |
| SO2 | 0.44 | 0.03 | 0.65 | Jun | Apr | 43 % |

(The components aren't independent, so the shares don't sum to 1.)

- **CO** has a very regular annual cycle: high from Dec to Apr (~40–44 mmol/m²) and low in the SW monsoon (~27 mmol/m² in Jul),
  repeating almost identically every year (`fig06`, `fig07`). It is the most predictable gas, mainly through its seasonal cycle.
- **NO2 and SO2** are dominated by window-to-window variability: the residual is about 60 % of the variance. Their seasonal cycles are
  weaker and less regular from year to year.

**Trend tests** use the seasonal Mann–Kendall test with the seasonal Sen slope on monthly area means (`trend_tests.csv`):

| Gas | 2019–2024 | p | Without 2020 | p |
|---|---|---|---|---|
| NO2 | +2.0 %/yr | 0.12 | +0.3 %/yr | 0.94 |
| CO | +0.4 %/yr | 0.36 | +0.6 %/yr | 0.29 |
| SO2 | −0.5 %/yr | 0.79 | −1.7 %/yr | 0.08 |

**No gas has a statistically significant trend over 2019–2024.** None of the 144 pixels has a significant trend
after FDR correction either (`fig11`). The apparent NO2 rise is mostly the 2020 dip followed by recovery:
dropping 2020 removes it. The annual means (`fig09`) show NO2 at −16 % in 2020 compared with 2019, then +14 % in 2023 and −3 % in 2024.
Six years is short for a trend test, so this means "no detectable trend", not "no trend".
SO2 shows a significant +5 %/yr over the train years alone (p = 0.014). That follows the same pattern of a 2020 low and then recovery,
so the slope depends on the period chosen and shouldn't be read as a long-term trend.

### COVID-19 lockdown (25 Mar – 31 May), a natural experiment (`fig12`)

| Gas | 2020 vs the same dates in other years |
|---|---|
| **NO2** | **−35.5 %** |
| SO2 | −21 % |
| CO | −0.7 % |

The NO2 drop is real: NO2 has no zeros in that period in any year, so missing data can't explain it. It shows that the
satellite NO2 signal responds to human activity (traffic and industry), which supports using it for forecasting.
CO doesn't respond: it is long-lived and dominated by the regional background and its seasonal cycle.

### Relationships between gases (`gas_correlation_spearman.csv`)

Spearman correlations of monthly anomalies (seasonal cycle removed): NO2–SO2 **0.51**, NO2–CO 0.23,
CO–SO2 ≈ 0. NO2 and SO2 vary together, consistent with shared combustion sources. CO behaves independently.

## 4. Spatial patterns (step 04_spatial)

- **NO2 has a clear urban–industrial hotspot** (`fig13`, `fig15`). The long-term maximum lies on the Manali–Ambattur axis in
  north-central Chennai, and pollution falls off strongly toward the south and the sea. The Manali/CPCL pixel ranks 4th of 144,
  at 1.68× the area median. Ambattur ranks 8th (1.58×) and Guindy 19th (1.42×) (`hub_enhancement.csv`).
- **The SO2 maximum is in the far north-east**, around 13.25–13.3°N, 80.35–80.4°E. That's north of the Ennore TPS hub point,
  near the coast where the Ennore/North-Chennai power stations and port are. Guindy, with lighter industry, is below the area median (0.91×).
- **CO is almost flat in space**: the long-term mean ranges only 36.2–37.6 mmol/m² (±2 %). It has no usable local
  hotspot signal at 5 km, and its variation is regional and seasonal.
- NO2 varies most (highest coefficient of variation) over the sea in the north-east. Plumes are carried offshore there,
  so concentrations swing between near-background and high.

## 5. Sentinel-2 (step 05_sentinel2)

- The median cloud-free share of an available S2 image is 75 % of the area.
- Area-median NDVI, NDBI and NDMI show no significant trend over 2019–2024 (`s2_index_trends.csv`, p ≥ 0.35). Land cover
  is effectively static at this scale over six years. S2 is better used as **static spatial context** (built-up density,
  vegetation, water) than as a time-varying input.
- The mean NDBI built-up pattern (`fig18`) matches the NO2 hotspot region qualitatively. The next step is to test this per S5P
  pixel.
- About 40 % of the raster is sea, so area-wide S2 statistics mix land and water. Later steps should mask the sea
  (NDVI < 0 together with NDMI > 0 is a reasonable first water mask, still to be checked).

## 6. What this means for forecasting

1. **Treat zeros as missing** and carry a validity mask. Never feed zeros to a model as low pollution.
2. **Remove the seasonal cycle.** It is the largest predictable component, especially for CO. Forecasting the anomaly against a
   per-pixel climatology fitted on the train years is the natural design, and the climatology itself is the baseline to beat.
3. **The 2020 lockdown is an anomaly inside the training years.** Consider flagging or down-weighting Mar–May 2020, or
   at least test whether excluding it changes the results.
4. **SO2 is noisy and 21 % of its windows are unusable**, so expect little skill beyond climatology. Report it, but don't let it drive design.
5. **NO2 is the most interesting target.** It has a strong, stable spatial hotspot, responds clearly to activity (the lockdown),
   has large short-term variability left to explain, and correlates with SO2.
6. **Sentinel-2 belongs as static context**, aggregated to the 12×12 S5P grid. It is missing in about half the monsoon windows and
   shows no change over time.

### Suggested next analyses
- Temporal autocorrelation of the anomalies, which tells us how far ahead persistence carries information and sets the forecast horizon.
- Per-pixel relationship between S2 land cover (built-up share, NDVI) and the long-term NO2 and SO2 means.
- Sensitivity of the SO2 results to the zeros-as-missing choice.
