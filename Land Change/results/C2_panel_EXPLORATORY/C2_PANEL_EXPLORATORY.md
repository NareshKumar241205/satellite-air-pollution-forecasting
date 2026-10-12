# C2 panel: EXPLORATORY (G1 failed at 0.5 %/pp; MDE 0.73 %/pp); raw and ERA5-deweathered Jan-Apr NO2

G1 detection limit (pooled long difference, 80 % power): **0.73 % NO2 per pp**. Estimates below that size can't be told apart from zero by this design; a null here is not evidence of no effect.
Sites: 18, years 2019-2023 (Jan-Apr). Spatial-block bootstrap: 500 draws of 4x4-pixel blocks within site.

| outcome                                            | model                                            | term                         |   estimate |   boot_se |   ci95_lo |   ci95_hi |   n_draws |
|:---------------------------------------------------|:-------------------------------------------------|:-----------------------------|-----------:|----------:|----------:|----------:|----------:|
| raw Jan-Apr NO2 (log of seasonal mean)             | panel two-way FE + SLX                           | L (own built-up)             |     -0.004 |     0.041 |    -0.081 |     0.085 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | panel two-way FE + SLX                           | WL (neighbour built-up)      |      0.267 |     0.108 |     0.049 |     0.476 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2020-2022 outcomes)                    | L                            |      0.039 |     0.055 |    -0.067 |     0.142 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2020-2022 outcomes)                    | L_lead (placebo: must be ~0) |      0.308 |     0.083 |     0.136 |     0.452 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2020-2022 outcomes)                    | L_lag                        |     -0.146 |     0.114 |    -0.354 |     0.099 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2020-2022 outcomes)                    | WL                           |      0.189 |     0.15  |    -0.093 |     0.453 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | long difference 2019->2023 + site FE (G1 design) | dL                           |     -0.026 |     0.064 |    -0.151 |     0.102 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | panel two-way FE + SLX                           | L (own built-up)             |     -0.002 |     0.041 |    -0.08  |     0.074 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | panel two-way FE + SLX                           | WL (neighbour built-up)      |      0.408 |     0.099 |     0.202 |     0.6   |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2020-2022 outcomes)                    | L                            |     -0.023 |     0.058 |    -0.13  |     0.095 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2020-2022 outcomes)                    | L_lead (placebo: must be ~0) |      0.35  |     0.069 |     0.187 |     0.459 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2020-2022 outcomes)                    | L_lag                        |      0.074 |     0.098 |    -0.141 |     0.249 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2020-2022 outcomes)                    | WL                           |      0.388 |     0.127 |     0.104 |     0.601 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | long difference 2019->2023 + site FE (G1 design) | dL                           |      0.003 |     0.068 |    -0.125 |     0.142 |       500 |

- raw Jan-Apr NO2 (log of seasonal mean): spatial placebo (500 draws) median b = 0.010, 95 % of |b| below 0.073; permutation p for the panel b = 0.922.
- ERA5-deweathered Jan-Apr NO2 (mean of window logs): spatial placebo (500 draws) median b = 0.011, 95 % of |b| below 0.074; permutation p for the panel b = 0.968.
