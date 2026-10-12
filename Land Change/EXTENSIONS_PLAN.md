# Extensions: rules frozen 13 Oct 2026, before any extension result

These are the optional steps of `PLAN.md` §4 that can run from this PC. Each is **exploratory or validation-only**. Nothing here re-opens frozen decisions:
- B = 0 for every gas;
- the G1 pooled long difference 2019→2023 stays the official gate;
- **2025 is spent** (each gas was scored once). No extension reads 2025 results for tuning, and none re-scores 2025.

A positive finding in any extension is a hypothesis for 2026 data, not a replacement result.

## X1. Lockdown-free land window (G1b)

- **Design.** Same data, 18 sites, Dynamic World label, same estimator (pooled long difference, pixel and site-year effects, spatial-block bootstrap). Window **2021→2024**, so 2020 is never an outcome or a land year. 2024 is validation in the forecast, not a test year.
- **Caveat.** The v1→v2 processor change (2021→2022) lies inside the window. City-year effects absorb a uniform step only.
- **Reported:** detectable effect at 0.5 and 1.0 %/pp; own effect with CI; lead test (outcomes 2022–2023).
- **Rule.** Exploratory. It can't change B. If G1b passes at 0.5 %/pp **and** its lead and neighbour terms are clean, the 2021→2025 check is run once as a confirmatory test and recorded whatever it shows. Otherwise it is not run.

## X2. CNN check (validation only)

- **Design.** The pre-registered small CNN (`forecast_run.cnn_fit`: 12×12 site grids, anomaly lags 0–5 plus missing flags, early stopping on validation). Run for NO2 and CO, 2019–2022 train, 2023–2024 validation.
- **Rule.** Compare validation MSE with the chosen Ridge. If the CNN wins, it is recorded as "validation winner, test deferred to 2026". The 2025 scores already recorded **stand**. No 2025 re-scoring. If PyTorch can't be installed, the step is reported as not run.

## X3. Interval bands: recalibration evidence (validation only)

- **Problem.** CO's 2025 coverage was 0.74 and SO2's 0.94 against 0.80.
- **Comparison, per gas and lead, on the persistence-based forecast Â(h) = φ_h·a_last with g = 1** (a documented simplification: the ML model only matters at lead 1):
  - **M0 (current):** global 10/90 % quantiles of the validation residual, fitted on 2023, evaluated on 2024.
  - **M1:** residual divided by the pixel's training-period anomaly sd; global quantiles of the standardized residual, fitted on 2023, bands = Â ± sd_p·q, evaluated on 2024.
- **Rule.** M1 is adopted for future forecasts only if its absolute coverage error on 2024 is smaller than M0's for at least two of the three gases at lead 1. Not applied retroactively to 2025.

## X4. Emission proxy: VIIRS night-time lights

- **Variable.** Annual VIIRS night-time lights radiance per block (`NOAA/VIIRS/DNB/ANNUAL_V22`, 2019–2024), as L = 10·ln(1 + radiance), so a 1-unit change is about a 10 % change in lights. No smoothing, no season (the product is annual).
- **Design.** Same estimator and window 2019→2023 (and 2021→2024), 18 sites. The G1 threshold is not defined for this unit, so only the detectable effect (% NO2 per 10 % change in lights) is reported, plus the own-effect estimate, lead and neighbour terms.
- **Rule.** Exploratory. Hypothesis only. It can't change B.

## X5. Climate indices for long-lead skill

- **Data.** ONI (ENSO) and DMI (Indian Ocean Dipole), monthly, NOAA public files.
- **Model.** For each lead h = 2…30, ridge regression of the block anomaly at t+h on [ONI(t), DMI(t)] with a fixed α, one site-specific coefficient pair per site, fitted on 2019–2022, evaluated on 2023–2024. Compared with climatology and with the validation-fitted damped persistence.
- **Rule.** Report the mean skill over leads 2–30 with a 95 % block-bootstrap CI. Adopt for the 2026 test only if the CI excludes 0 and the skill is > 0. Otherwise: "no long-lead skill from climate indices".

## Not executable from this PC

- Consistent reprocessed NO2 (v2.4.0 Level-2): needs a Copernicus account and ~300 GB.
- Pull request into `main`: `gh` isn't logged in.
- Deck cover placeholders and chapter framing: need Naresh.
