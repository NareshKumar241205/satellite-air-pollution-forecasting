# Site selection rule for the added G1 sites (frozen 12 Oct 2026, before any NO2–land slope is estimated)

Decision (Naresh, via the project thread): G1 is the **pooled long difference 2019→2023**, powered for
**0.5 % NO2 per pp** (also reported at 1.0 %/pp). Sites are added **only by Dynamic World built-up growth, never
by NO2**. This file is committed before the selection runs and before any slope is looked at. The rule must not be
changed after seeing G1 or any slope.

## What has been seen so far (for the record)

- **Looked at:** the G1 noise statistics (sd of the city-year-demeaned dNO2 and its neighbour correlation), and the land-change distributions of the 18 current cities.
- **Not looked at:** no regression of NO2 on land has been run, and no NO2–land correlation has been computed.

## Rule

1. **Candidate tiles.** Tile the TNAQ export grid (0.01°, origin 76.22 E / 13.54 N) into 48 × 48 sub-pixel tiles
   (12 × 12 blocks of 0.04°). Tile corners are at multiples of 48 sub-pixels from the origin; partial tiles at the
   edges are dropped.
2. **Geometry filters (no NO2 values):**
   - drop tiles that overlap any existing site box (`configs/land.yaml` cities);
   - keep a tile only if ≥ 50 % of its blocks are inside the export footprint.

   The footprint is the set of sub-pixels valid in at least one 2019 window. It is the Tamil Nadu clip of the export, so it depends on data availability, not NO2 values.
3. **Land data.** Same-season (Jan–Apr) Dynamic World V1 mode-label fractions for 2019–2025, from
   `scripts/ee_s2_landfrac.py --classifier dynamicworld`, exactly as for the cities.
4. **Growth score.** The median over the tile's land blocks of built-up(2023) − built-up(2019), in pp. Land blocks are those with DW `valid_frac` ≥ 0.5 and water fraction < 0.5 in 2019.
5. **Selection.** Add every tile with a growth score ≥ **1.5 pp**. The threshold is fixed a priori: it is just above the
   pooled median long-difference gain of the 18 current cities (1.39 pp), so added sites grow at least as fast as the current pool.
   All qualifying tiles are added; there is no cap and no ranking.
6. The selected tiles are written to `configs/sites_added.yaml` (name `T<row>_<col>` plus centre lat/lon) and
   are loaded as extra sites by `aqf_land.tnaq.load_land_config`.

Nothing in steps 1–5 reads NO2 values. Implemented in `src/aqf_land/sites.py`.
