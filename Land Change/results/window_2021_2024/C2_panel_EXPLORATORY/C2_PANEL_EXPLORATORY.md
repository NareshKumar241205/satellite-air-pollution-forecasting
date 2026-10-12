# C2 panel: EXPLORATORY (G1 failed at 0.5 %/pp; MDE 0.73 %/pp); raw and ERA5-deweathered Jan-Apr NO2

G1 detection limit (pooled long difference, 80 % power): **1.29 % NO2 per pp**. Estimates below that size can't be told apart from zero by this design; a null here is not evidence of no effect.
Sites: 18, years 2021-2024 (Jan-Apr). Spatial-block bootstrap: 500 draws of 4x4-pixel blocks within site.

| outcome                                            | model                                            | term                         |   estimate |   boot_se |   ci95_lo |   ci95_hi |   n_draws |
|:---------------------------------------------------|:-------------------------------------------------|:-----------------------------|-----------:|----------:|----------:|----------:|----------:|
| raw Jan-Apr NO2 (log of seasonal mean)             | panel two-way FE + SLX                           | L (own built-up)             |     -0.056 |     0.042 |    -0.132 |     0.033 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | panel two-way FE + SLX                           | WL (neighbour built-up)      |     -0.24  |     0.126 |    -0.486 |     0.01  |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2022-2023 outcomes)                    | L                            |     -0.083 |     0.121 |    -0.265 |     0.201 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2022-2023 outcomes)                    | L_lead (placebo: must be ~0) |      0.157 |     0.135 |    -0.094 |     0.439 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2022-2023 outcomes)                    | L_lag                        |     -0.346 |     0.15  |    -0.549 |     0.024 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2022-2023 outcomes)                    | WL                           |      0.145 |     0.238 |    -0.317 |     0.625 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | long difference 2021->2024 + site FE (G1 design) | dL                           |     -0.188 |     0.07  |    -0.318 |    -0.045 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | panel two-way FE + SLX                           | L (own built-up)             |     -0.071 |     0.047 |    -0.157 |     0.022 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | panel two-way FE + SLX                           | WL (neighbour built-up)      |     -0.113 |     0.094 |    -0.287 |     0.073 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2022-2023 outcomes)                    | L                            |      0.112 |     0.107 |    -0.092 |     0.317 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2022-2023 outcomes)                    | L_lead (placebo: must be ~0) |      0.452 |     0.147 |     0.137 |     0.706 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2022-2023 outcomes)                    | L_lag                        |     -0.053 |     0.137 |    -0.308 |     0.206 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2022-2023 outcomes)                    | WL                           |      0.153 |     0.18  |    -0.228 |     0.498 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | long difference 2021->2024 + site FE (G1 design) | dL                           |     -0.105 |     0.081 |    -0.256 |     0.057 |       500 |

- raw Jan-Apr NO2 (log of seasonal mean): spatial placebo (500 draws) median b = 0.005, 95 % of |b| below 0.105; permutation p for the panel b = 0.288.
- ERA5-deweathered Jan-Apr NO2 (mean of window logs): spatial placebo (500 draws) median b = 0.013, 95 % of |b| below 0.104; permutation p for the panel b = 0.178.
