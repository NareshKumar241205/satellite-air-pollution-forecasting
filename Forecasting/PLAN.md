# Forecasting plan

Two forecasting tasks run on the same data, the same rules and the same evaluation:

| Task | Lead | Output |
|---|---|---|
| **A. Short-term** | t+1: the next 5-day window | the 12×12 map per gas, with 10/90 % bands |
| **B. Long-term outlook** | t+1 … t+30: 5 to 150 days | 30 maps per gas, 10/90 % bands, area-mean outlook |

The rules follow from `Data Analysis/REPORT.md`:
- Exact zeros are missing. They are masked out of every loss and every metric, never used as a target.
- Split by year: **train 2019–2022, validate 2023, test 2024**. A (origin, lead) pair belongs to the split of its
  *target* window. The climatology, scalers and every fitted parameter come from train years only.
- Models predict the **anomaly** against a per-pixel harmonic climatology fitted on train. Metrics are reported in mol/m².
- 2024 is scored **once**, at the end (F5), with frozen configs. Only after that is anything refit on 2019–2024
  for the operational outlook (Jan–May 2025).

## Phases

| Phase | Content | Status |
|---|---|---|
| **F0 Forecastability** | anomaly autocorrelation for leads 1–30 (pooled pixels and area mean); cross-gas lead–lag partial correlation; neighbour (3×3) predictive correlation; Sentinel-2 NDBI vs long-term NO2 per S5P pixel with a spatial-block permutation test | done, see REPORT.md |
| **F1 Baselines** | climatology, persistence, damped anomaly persistence; task A and task B; RMSE, MAE, skill vs climatology, anomaly correlation; residual-quantile 10/90 % intervals with coverage and CRPS per lead; block-bootstrap CIs | done, see REPORT.md |
| F2 Learned models | per-gas LightGBM pooled over pixels plus a Ridge reference. Inputs are anomaly lags, the 3×3 neighbourhood, other gases, season, static S2 per pixel and validity masks. For task B, the lead is a feature (direct multi-horizon). Quantile versions give the 10/90 % bands | done, see REPORT.md |
| F3 ConvLSTM | small spatio-temporal network, **only if** F0/F2 show spatial/temporal structure that the tree models miss | not justified (F2), skipped |
| F4 Ablations | S2 off / S2 shuffled; own gas vs multi-gas; no neighbours; no validity masks; exclude 2020 Q2; plus rolling-origin years 2021–2023 and an NE-monsoon gate | done, see REPORT.md (model choice revised) |
| F5 Test | single scoring of the 2024 test year with frozen configs, then refit on 2019–2024 and produce the Jan–May 2025 outlook | waiting for go-ahead |

## Deliverables (F5)
- Task A: a GeoTIFF of the t+1 map per gas with its 10/90 % bands, and a CSV of the area mean.
- Task B: a 30-band GeoTIFF per gas (mean, 10 %, 90 %), a 30-row CSV outlook, hotspot series (Manali NO2, Ennore SO2) and fan charts.
- A skill-vs-lead curve per gas, with interval coverage and CRPS, on validation and on test.

Results go to `Forecasting/results/<phase>/`. The written findings go to `Forecasting/REPORT.md`.
