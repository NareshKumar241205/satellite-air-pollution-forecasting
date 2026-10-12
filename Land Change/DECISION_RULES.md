# Decision rules and frozen forecasting settings (committed 12 Oct 2026, before the results they govern)

Naresh's standing instruction: plan and execute end to end, then report once. Every rule below is fixed before its
result is seen. Nothing here may be changed after looking at validation or test numbers, except where a rule
itself says how a choice is made on validation.

## B. Land term in the forecast

**Steadier layer.** It is chosen on noise and lead-test criteria only (`DIAGNOSTICS_PLAN.md`):
- The candidates are `dw_prob_trail3` and `obt`.
- Pick the one with the smaller G1 MDE among those whose lead test is clean.
- If neither is clean, there is no steadier layer and B = 0.

**Rule.** If G1 passes at 0.5 %/pp (pooled long difference 2019→2023) **and** the lead test is clean on the steadier layer:
- run the confirmatory C2 panel once on that layer;
- set B[p, y] = β̂·L[p, y] + θ̂·WL[p, y] (log-space, as a multiplicative adjustment relative to the training-period mean land).

Otherwise **B = 0**. The thesis result is then the detection limit: "no land→NO2 effect above X %/pp is detectable".

"Clean" means two things:
- the lead coefficient's 95 % bootstrap CI covers 0 on all outcome years 2020–2022, with no years dropped;
- the neighbour term's CI covers 0.

## C. Forecasting (TNAQ sites, 0.04° blocks, 5-day windows, NO2 only)

**Data and splits.**
- **Target:** the block-mean tropospheric NO2 of the window, as an anomaly from the per-pixel harmonic climatology (K = 3, fitted on train).
- **Splits:** train 2019–2022 (windows starting in those years), validation 2023–2024, test 2025. 2025 is used once, in step E.

**Short range (C3, lead h = 1, the next window).** This is a ladder of models.
1. **Climatology.**
2. **Damped persistence.** a(t+h) = φ_h·a_last, with φ_h fitted on train and a_last the last valid anomaly within 3 windows (otherwise 0).
3. **Ridge** (α picked on validation from {0.1, 1, 10, 100, 1000}).
4. **LightGBM.** Settings are fixed: learning_rate 0.05, num_leaves 31, min_data_in_leaf 200, feature_fraction 0.8, bagging 0.8/1, λ2 1, up to 2000 rounds with early stopping 50 on validation.

**Features (all known at the origin).**
- own anomaly lags 0–5 and missing flags;
- 3×3 neighbour mean anomaly, lags 0–2;
- site mean anomaly, lags 0–2;
- trailing means of 6, 12 and 18 windows;
- day-of-year harmonics;
- ERA5 site scalars of windows t and t−1, standardised on train (past only);
- land: built-up fraction from the most recent complete Jan–Apr season strictly before the origin window (year y−1 for origins before May of year y, year y otherwise).

**Optional CNN.** A small CNN or ConvGRU is tried on 12×12 site grids only if PyTorch installs. It is used only if it beats LightGBM on validation MSE.

**Model choice for h = 1.** The lowest validation MSE among the ladder (and the CNN if run) wins.

**Long range (h = 1..30).** Climatology and damped persistence only, giving a skill-vs-lead curve.

**Metrics.**
- skill vs climatology, 1 − MSE/MSE_clim;
- RMSE, MAE and pattern correlation (the spatial correlation of anomalies per site-window, averaged);
- per-season skill (IMD seasons);
- 95 % block-bootstrap CIs over origins, in blocks of 6 consecutive windows;
- coverage of the 10–90 % band, taken from validation residual quantiles per lead.

**Ablations (validation, LightGBM).** ERA5 off; land off; land shuffled across pixels within site.

## D. Fusion

- F(h) = Clim + B + g(h)·Â(h).
- Â(h) is the chosen h = 1 model at lead 1 and damped persistence at h > 1.
- g(h) is fitted per lead on validation by least squares, clipped to [0, 1], then made monotone non-increasing in h (isotonic).

## E. Final scoring

- Freeze everything above.
- Refit the climatology, φ_h, the chosen model (with its validation-chosen settings and round count), g(h) and the interval quantiles on train + validation (2019–2024).
- Score 2025 once, and write down the numbers whatever they are.

## G. Other gases: CO and SO2 (added 13 Oct 2026, committed before any CO or SO2 result)

Scope added at Naresh's request ("update these tasks too"). **Every rule in sections C–E applies unchanged, gas by gas.**

**Same.**
- sites, blocks, windows, splits and seed;
- the features, models and fixed settings;
- the metrics, the bootstrap and the fusion rule;
- land term B = 0 for both gases. The land decision for NO2 was made on NO2 and is not re-opened. The land power check for CO and SO2 is run (section G.3) only as a labelled exploratory item.

**Differences.**
1. **Data.**
   - **Source:** the TNAQ export for 2019–2024, plus 2025 downloaded with the same repaired method (`ee_s5p_no2_windows.py --gas`). Bit-for-bit reproduction of two existing windows per gas is the acceptance check, and the NaN, −9999 and valid classes must also match. If the check fails, 2025 is not scored for that gas.
   - **Units:** µmol/m² for SO2, mmol/m² for CO. Zeros are missing, as in the main data rules (`configs/data.yaml`).
2. **ERA5 and land features.** They are the same, so each gas gets the same feature set.
3. **Land power check, exploratory.** Gate G1 and the C2 panel are re-run per gas on the DW-label layer only, 2019–2023 tuning years, with the unchanged estimator. Results are labelled exploratory and cannot change any decision.
4. **Test.** Each gas is scored on 2025 exactly once after its own validation stage, and the number is recorded whatever it is. The NO2 test is not re-scored.

**Reading rule.** If a gas's validation skill at lead 1 has a 95 % CI that covers 0, the gas is reported as "no skill beyond climatology" and the climatology is its forecast. This is stated before any result is known.
