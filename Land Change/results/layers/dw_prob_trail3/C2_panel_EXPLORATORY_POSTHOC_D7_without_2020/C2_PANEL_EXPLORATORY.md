# C2 panel: POST-HOC DIAGNOSTIC D7, outcome years [2020] removed (DIAGNOSTICS_PLAN.md addendum). EXPLORATORY (G1 failed at 0.5 %/pp; MDE 0.73 %/pp); raw and ERA5-deweathered Jan-Apr NO2

G1 detection limit (pooled long difference, 80 % power): **1.52 % NO2 per pp**. Estimates below that size can't be told apart from zero by this design; a null here is not evidence of no effect.
Sites: 18, years 2019-2023 (Jan-Apr). Spatial-block bootstrap: 500 draws of 4x4-pixel blocks within site.

| outcome                                            | model                                            | term                         |   estimate |   boot_se |   ci95_lo |   ci95_hi |   n_draws |
|:---------------------------------------------------|:-------------------------------------------------|:-----------------------------|-----------:|----------:|----------:|----------:|----------:|
| raw Jan-Apr NO2 (log of seasonal mean)             | panel two-way FE + SLX                           | L (own built-up)             |     -0.037 |     0.08  |    -0.191 |     0.121 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | panel two-way FE + SLX                           | WL (neighbour built-up)      |      0.231 |     0.236 |    -0.261 |     0.699 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2021-2022 outcomes)                    | L                            |     -0.219 |     0.203 |    -0.652 |     0.177 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2021-2022 outcomes)                    | L_lead (placebo: must be ~0) |      0.247 |     0.389 |    -0.486 |     1.03  |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2021-2022 outcomes)                    | L_lag                        |     -0.373 |     0.322 |    -0.915 |     0.281 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2021-2022 outcomes)                    | WL                           |     -0.085 |     0.506 |    -1.098 |     0.895 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | long difference 2019->2023 + site FE (G1 design) | dL                           |     -0.065 |     0.101 |    -0.275 |     0.126 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | panel two-way FE + SLX                           | L (own built-up)             |     -0.036 |     0.075 |    -0.187 |     0.098 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | panel two-way FE + SLX                           | WL (neighbour built-up)      |      0.576 |     0.231 |     0.109 |     0.989 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2021-2022 outcomes)                    | L                            |     -0.461 |     0.177 |    -0.858 |    -0.136 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2021-2022 outcomes)                    | L_lead (placebo: must be ~0) |      0.102 |     0.352 |    -0.554 |     0.81  |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2021-2022 outcomes)                    | L_lag                        |     -0.274 |     0.32  |    -0.827 |     0.417 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2021-2022 outcomes)                    | WL                           |      0.075 |     0.494 |    -0.894 |     1.075 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | long difference 2019->2023 + site FE (G1 design) | dL                           |      0     |     0.107 |    -0.217 |     0.192 |       500 |

- raw Jan-Apr NO2 (log of seasonal mean): spatial placebo (500 draws) median b = -0.012, 95 % of |b| below 0.158; permutation p for the panel b = 0.664.
- ERA5-deweathered Jan-Apr NO2 (mean of window logs): spatial placebo (500 draws) median b = -0.002, 95 % of |b| below 0.156; permutation p for the panel b = 0.640.
