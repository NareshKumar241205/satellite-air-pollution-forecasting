# Chapter X: Land-use change, Sentinel-5P NO2 and short- to long-range forecasting over Tamil Nadu

*Draft for the thesis. All numbers are from `Land Change/results/` and are summarised in `docs/FINAL_REPORT.md`.*

## X.1 Introduction

Urban expansion is often assumed to raise nitrogen dioxide (NO2), the main tracer of combustion from traffic and
industry. Satellite instruments such as TROPOMI on Sentinel-5P measure the tropospheric NO2 column daily at about
5 km, and Sentinel-2 maps land cover at 10 m. Together they make it possible to ask, for every 5 km pixel and every
year, whether places that gained built-up area also gained NO2. A second question is practical: how far ahead can
NO2 maps be forecast, and does land-cover information improve the forecast?

The two clocks of the problem shape the design. Land changes slowly, by a few percentage points over several
years. NO2 changes from day to day with emissions, weather and chemistry. This chapter measures land once per year
from a same-season (January–April) Sentinel-2 composite. It compares that with the same-season NO2 of the same year,
and it treats weather and retrieval changes as nuisance factors common to a city and year.

## X.2 Data and study design

The study uses 18 non-overlapping 0.48° boxes centred on Tamil Nadu cities and towns. Each box is divided into
12 × 12 blocks of 0.04° (about 4.4 km), close to the TROPOMI footprint, giving 2,568 blocks.

- **NO2:** five-day composites of the operational TROPOMI tropospheric column (Earth Engine `COPERNICUS/S5P/OFFL/L3_NO2`), 2019–2025.
- **Land:** built-up fractions per block and year from Dynamic World (10 m, Sentinel-2 based). Two smoother alternatives were tested: a trailing three-year median of the Dynamic World built-up probability, and Open Buildings Temporal.
- **Weather:** ERA5 hourly fields around the satellite overpass, summarised per city and five-day window.
- **Retrieval versions:** the TROPOMI processor version was checked image by image. Every January–April season from 2019 to 2025 used a different version, and the largest change (v1 to v2) fell between 2021 and 2022.

The decision rules were written and committed before the results they govern:
- the site-selection rule;
- the diagnostic plan;
- the decision on whether land enters the forecast;
- every forecasting setting.

Years 2019–2023 were used for the land analysis. For forecasting, 2019–2022 were used for training and 2023–2024 for validation, and 2025 was held back and scored only at the end.

## X.3 Can land-use change be detected in NO2?

### X.3.1 Power before estimation

Before any land–NO2 slope was estimated, a power check (gate G1) asked what effect the design could detect. G1 combined three measured quantities:
- the year-to-year noise in NO2 after removing city-year effects;
- the spatial correlation of that noise between neighbouring blocks;
- the real spread of land change after removing land-label noise.

On the best design, the 2019→2023 long difference over all 18 sites with the Dynamic World label, the minimum detectable effect was **0.73 % NO2 per percentage point of built-up gain**. That is above the 0.5 %/pp the study set out to detect, so G1 failed. The design could only detect effects of 1 %/pp or more.

Two attempts to improve power did not help:
- **More sites.** Adding sites selected only for fast built-up growth found no further qualifying area in the state.
- **Smoother layers.** These reduced label noise, but they reduced the measured change even more, and the detectable effect rose to 1.5 and 4.2 %/pp.

An early check on a 1 km Sentinel-2 index proxy was also instructive. Its "built-up change" correlated −0.99 with vegetation change: it was tracking rainfall and the 2019 drought, not urbanisation.

### X.3.2 Exploratory estimates and why they cannot be read causally

**Model.** A panel model regressed log NO2 on built-up fraction, with block and city-year fixed effects and a neighbour term. Uncertainty came from a spatial-block bootstrap.

**Estimates.** The own-block estimate was close to zero on every layer and outcome: −0.004 %/pp (95 % CI −0.08 to 0.09) with the Dynamic World label. The long difference agreed.

**The lead test failed.** Next year's built-up fraction "predicted" this year's NO2 (+0.31 %/pp, CI 0.14–0.45), which a causal effect cannot do. The neighbour term was also significant.

**Diagnostics.** A pre-registered set of diagnostics traced the failure:
- **Lockdown year.** Removing the 2020 outcome year eliminated it (0.31 → 0.01). January–April 2020 contains the COVID-19 lockdown, whose NO2 reduction varied across space; the most plausible reading is that this variation was related to later growth, although the diagnostic shows only that 2020 carries the lead.
- **Aerosol.** In 2019–2021, the years for which the aerosol data could be retrieved, aerosol optical depth correlated strongly with the built-up label (r = 0.61), and next-year aerosol correlated with both the land lead and the NO2 residual.
- **Ruled out:** single sites, box edges and shared spatial residual structure.

**Conclusion of this section.** At ~5 km over 2019–2023, the satellite data do not show a land→NO2 effect larger than about 0.7 % per percentage point. The available contrasts are also confounded by the lockdown year and by aerosol. Under the pre-registered rule, land therefore enters the forecast with zero weight.

## X.4 Forecasting NO2

The forecast target is the NO2 anomaly of each block relative to its own seasonal climatology (three annual harmonics, fitted on the training years). The models form a ladder:
1. climatology;
2. damped persistence;
3. ridge regression;
4. gradient-boosted trees.

Features were known at forecast time:
- recent NO2 anomalies of the block, its neighbours and its city;
- trailing means;
- the time of year;
- ERA5 weather;
- the latest built-up fraction.

The model was chosen on validation. A small convolutional network was planned but not run, because the deep-learning library could not be installed.

**Validation (2023–2024).** For the next five-day window, ridge regression improved the mean squared error by **9.4 %** over climatology (CI 6.4–13.4 %). Gradient boosting gained 6.5 % and damped persistence 5.4 %. Weather added a little skill. Removing or shuffling land fractions changed nothing, consistent with Section X.3.

**Test (2025).** With every setting frozen and refitted on 2019–2024, the five-day forecast improved on climatology by **3.6 %** (CI −1.3 to 6.8 %). The interval includes zero, so the test-year gain is positive but not statistically significant. The 10–90 % band contained 81 % of observations. From ten days to 150 days, skill was indistinguishable from zero. The skill-versus-lead curve drops from the first window to the second and stays flat.

**Long range.** The useful long-range forecast is the seasonal climatology itself, with calibrated uncertainty. At the city level, 78 % of January–May 2025 site means fell inside the 10–90 % band, and those site means ran a median 7 % below the 2019–2024 climatology.

## X.5 Discussion

The negative land result has three likely causes:
- **Signal versus noise:** a 4 km block changes by only 1–2 percentage points of built-up area over four years, comparable to the noise in global land-cover labels.
- **Retrieval changes:** TROPOMI retrievals changed version every year.
- **Lockdown confounding:** the lockdown happened in the middle of the land window.

A detection limit is still informative. It bounds how large an urbanisation effect could be at this scale. It also tells future work what is needed:
- a consistently reprocessed NO2 record;
- a longer land window that avoids 2020;
- higher-resolution or emission-weighted land variables (traffic, industry) instead of built-up area alone.

For forecasting, the result matches the earlier Chennai-only work: beyond a few days, NO2 at these scales is predictable mainly through its seasonal cycle.

## X.6 Integrity notes

The first 2025 scoring had to be discarded. Pixels outside Tamil Nadu had been filled with zeros instead of
missing values, which faked large NO2 drops at coastal and border sites. The fault was found from those drops and traced and fixed. 2025 was then scored again with unchanged settings, and both records are kept. Earlier in the work, a spurious power "pass" caused by label noise was found and corrected before any slope was estimated.
