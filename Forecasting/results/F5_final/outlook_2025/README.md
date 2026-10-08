# Jan–May 2025 outlook (F5)

Issued from the last observed window (25 Dec 2024). 30 steps of 5 days: 2024-12-30 to 2025-05-24. Everything was refitted on 2019–2024, with the climatology recomputed on those years.

Model used per lead (frozen before the 2024 test, never changed afterwards):

- **NO2**: leads 1-1 (5-5 d): damped; leads 2-6 (10-30 d): ridgeB; leads 7-30 (35-150 d): climatology
- **CO**: leads 1-1 (5-5 d): damped; leads 2-30 (10-150 d): climatology
- **SO2**: leads 1-30 (5-150 d): climatology

Files:
- `outlook_<GAS>_2025.tif`: bands 1–30 are the mean, 31–60 the 10 % quantile, 61–90 the 90 % quantile, in mol/m².
  Each band's description and tags give the lead, the window start and the model used.
- `outlook_area_mean.csv`: area-mean forecast with an 80 % band, one row per gas and lead.
- `outlook_hotspots.csv`: forecast at the Manali (NO2) and Ennore (SO2) hub pixels.
- `fig_outlook_fan_area.png`, `fig_outlook_hotspots.png`, `fig_outlook_maps.png`.

Interpretation: beyond about 30 days the forecast is the seasonal climatology (Forecasting/REPORT.md, F4).
The bands are 80 % intervals from leave-one-year-out residuals. Pixel bands are clipped at 0.
