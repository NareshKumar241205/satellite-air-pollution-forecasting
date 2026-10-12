# Lead-test diagnostics and steadier land layers (plan frozen 12 Oct 2026, before running any of it)

The exploratory C2 panel failed its lead test: next year's built-up L(t+1) "predicts" this year's NO2 (+0.31 %/pp
raw, +0.35 deweathered; `results/C2_panel_EXPLORATORY/`). This file fixes, in advance, what is checked and how the
outcome is read. Estimators, the site set, the 1.5 pp site rule, tuning years 2019–2023 and the reserved 2025 test
year don't change. No result below may be used to pick a different estimator.

## Part 1: diagnostics on the current layer (Dynamic World mode label)

Each check reports its numbers whatever they show.

| ID | Question | How | Read as an artifact if |
|---|---|---|---|
| D1 | Do DW labels move with scene or atmospheric conditions at the pixel level? | Per block and year (Jan–Apr): DW image count, mean S2 cloud probability (Cloud Score+ 1 − cs_cdf), mean MODIS MAIAC AOD 0.55 µm (`MODIS/061/MCD19A2_GRANULES`), and median S2 brightness (mean of B2–B4 reflectance). Within-site, within-year correlation of each with L(t), with L(t+1) − L(t), and with the NO2 residual. | A covariate correlates with both L(t+1) and NO2(t) at \|r\| ≥ 0.1 |
| D2 | Is the year-to-year label change level-dependent in particular years? | For each consecutive pair, the within-site slope of dL on the pixel's NO2 level and on its L level. | The slope differs sharply in one pair (e.g. one with an S2 or DW processing change) |
| D3 | Do pixel-specific trends explain the lead? | Re-fit the lead/lag panel with pixel-specific linear trends. | The lead coefficient drops to ≈ 0 (its 95 % CI covers 0) |
| D4 | Is it driven by one site or one year? | Leave-one-site-out and leave-one-outcome-year-out lead coefficients. | One site or year carries most of it (dropping it moves the lead by > 50 %) |
| D5 | Edge or spillover effects at box borders? | Drop the outer ring of blocks in every site and re-fit. | The lead changes by > 50 % |
| D6 | Is there a common spatial component? | Moran-type neighbour correlation of the NO2 residual and of the L(t+1) residual from the lead/lag model, and their cross-correlation between neighbours. | Strong shared neighbour structure (cross r ≥ 0.2) |

## Part 2: steadier built-up layers (defined here, before computing)

- **S1: DW built-up probability, trailing 3-year median.** Per block and year y, take the mean Jan–Apr DW `built` probability (10 m, mean over images, then over the block). The layer value is the median of years y−2, y−1, y. The window is **trailing, not centred**: a centred window would put L(t+1) into L(t) and create a lead correlation mechanically. Years 2017–2024 are fetched, so 2019 has its 3-year window.
- **S2: Open Buildings Temporal** (`GOOGLE/Research/open-buildings-temporal/v1`), if it covers the sites. Per block and year, use the mean `building_fractional_count`-based building presence. The annual inference is used as is; the product is already temporally smoothed by design. Used for the years it covers within 2019–2023.

Each layer replaces the built-up fraction only; NO2, sites, blocks and years stay the same. For each layer:
- G1: pooled long difference 2019→2023, MDE at 0.5 and 1.0 %/pp.
- The exploratory C2 panel: same model, same bootstrap and placebo, with the lead test.

Reported for each layer:
- the MDE;
- the lead coefficient with its CI;
- whether the confounding signature is gone (the lead CI covers 0 **and** the neighbour term's CI covers 0).

Whichever layer is used later, the choice must be made on the noise and lead-test criteria above, not on the size of the own-pixel slope.

## Addendum (13 Oct 2026, committed before running): post-hoc diagnostic D7

D4 showed that the lead is carried by outcome year 2020 (the lockdown). D7 re-runs the same C2 panel, with the same layers, bootstrap and placebo, with the **2020 outcomes removed**. 2020 land values are kept, so they can still act as leads or lags of other years. The long difference 2019→2023 is unchanged.

D7 is **post hoc and labelled as such**. It cannot change rule B (B = 0 stands, because G1 fails), the site set, or any forecast setting. It is reported to show what the panel says once the lockdown year is out:
- the lead coefficient (outcomes 2021–2022);
- the neighbour term;
- the own-pixel slope.
