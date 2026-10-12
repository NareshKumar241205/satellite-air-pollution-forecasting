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

### Processor versions (`results/P1_provenance/PROVENANCE.md`)

- **The local export doesn't record the TROPOMI processor version.** The GeoTIFFs carry only `AREA_OR_POINT`, and the manifests have no version column.
- By the public OFFL timeline (approximate), the Jan–Apr seasons of 2019–2021 are v1.x, and those of 2022–2024 are v2.3.1 or later. So the 2021→2022 pair straddles the v1→v2 change.
- The data shows no step at that pair:
  - The mean change across cities is −1.6 %. Other pairs range from −22 % (2020, lockdown) to +16 %.
  - The level dependence in 2021→2022 (−11 % per log-unit) is the same regression to the mean every pair shows (−2 to −10).
- So the step is either small over Tamil Nadu or absent from this collection. **Exact versions need `scripts/ee_s5p_versions.py`.**

### Gate G1 (`results/G1_power/`): **FAIL on the local data, and the cause is the land proxy, not the design**

| Pooled, 9 cities | pairs | sd(dL) | noise sd | signal sd | corr(dL, dVeg) | sd(dNO2 resid) | N_eff | needs true sd(dL) ≥ |
|---|---|---|---|---|---|---|---|---|
| consecutive years | 4,456 | 18.5 pp | 32.3 | **0** | **−0.99** | 7.8 % | 539 | **1.9 pp/yr** |
| without 2020 | 2,228 | 15.1 | 21.7 | 0 | | 7.4 % | 323 | 2.3 pp/yr |
| without 2021→22 | 3,342 | 19.2 | 33.6 | 0 | | 8.3 % | 363 | 2.4 pp/yr |
| one 2019→2023 difference | 1,114 | 24.3 | 46.7 | 0 | | 6.6 % | 155 | **3.0 pp over 4 yr** |

**The 1 km "built-up" dL is greenness, not land use.**
- Within city-year it correlates −0.99 with the vegetation change.
- Its negative side (built-up loss, which is rare in reality) is larger than the whole spread, so no signal remains after removing classification noise.
- From 2019 to 2023, "built-up" falls by more than 30 pp almost everywhere, because the 2019 dry season was a drought year (`fig_g1_cumulative_dL_maps.png`).

A first version of G1 that used the raw sd(dL) reported PASS (MDE 0.05 %/pp). That result was spurious, and it's why G1 now subtracts the noise.

**What is measured, and holds independently of the proxy, is the NO2 side:**
- the residual dNO2 noise after city-year fixed effects (5–16 % per city);
- its spatial correlation between neighbouring blocks (ρ 0.6–0.9).

From these, the last column gives the true built-up spread needed to detect 0.5 %/pp at 80 % power.

**Chennai alone can't work.** It needs a true sd(dL) ≥ 15 pp per year, or ≥ 15.5 pp over 4 years, because its NO2 residual is noisy and strongly autocorrelated (N_eff 23). Pooled, the 9 cities need ≈ 2 pp per year or ≈ 3 pp over 4 years. Peri-urban growth of a few pp over 4 years per 4 km block is plausible on the urban fringe but not across whole boxes. So the long difference, 9+ cities and a phenology-robust built-up layer together might pass, and the real 20 m / 10 m dL decides. **Recommendation:** keep all nine cities and add more if the EE dL spread is below the threshold; use the long difference rather than year-to-year pairs; use Dynamic World (`--classifier dynamicworld`) rather than an NDBI rule.

## What still needs credentials (not run, nothing was downloaded)

| Step | Script | Needs |
|---|---|---|
| Exact processor versions | `uv run --with earthengine-api python "Land Change/scripts/ee_s5p_versions.py" --project <gcp>` | Earth Engine account |
| Real dL: same-season S2 SR 20 m (or Dynamic World 10 m) land fractions per block and year | `uv run --with earthengine-api python "Land Change/scripts/ee_s2_landfrac.py" --project <gcp> --years 2019 2020 2021 2022 2023 2024 2025 --classifier dynamicworld`, then download the CSVs into `Land Change/data/ee_landfrac/` and re-run `uv run aqf-land all` (it switches to them automatically) | Earth Engine + Google Drive |
| ERA5 per year, padded box, overpass hours → scalars per city and window | `uv run --with cdsapi python "Land Change/scripts/cds_era5.py" download`, then `uv run --with xarray --with netcdf4 python "Land Change/scripts/cds_era5.py" reduce` | CDS account, `~/.cdsapirc`, ERA5 licence accepted |
| 2025 test year (NO2 + S2) | `AirQ/Scripts/Data_acquisition.js` with `END_DATE = '2025-12-31'`, `YEAR_FILTER = 2025` | Earth Engine Code Editor |

The ERA5 reducer was tested end-to-end on a synthetic NetCDF only; no ERA5 data exists on this machine.

## Files

- `src/aqf_land/tnaq.py`: city boxes, cached sub-pixel reads, block means.
- `src/aqf_land/layers.py`: seasonal composite, classification, class fractions, seasonal NO2, the pixel-year pair table, and the EE CSV loader.
- `src/aqf_land/provenance.py`, `g1.py`, `era5.py`, `cli.py`.
- `scripts/`: Earth Engine and CDS jobs.
- `results/P1_layers/` is gitignored (rebuildable); `P1_provenance/` and `G1_power/` are tracked.
