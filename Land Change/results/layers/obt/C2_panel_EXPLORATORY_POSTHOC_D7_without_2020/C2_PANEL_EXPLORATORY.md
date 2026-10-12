# C2 panel: POST-HOC DIAGNOSTIC D7, outcome years [2020] removed (DIAGNOSTICS_PLAN.md addendum). EXPLORATORY (G1 failed at 0.5 %/pp; MDE 0.73 %/pp); raw and ERA5-deweathered Jan-Apr NO2

G1 detection limit (pooled long difference, 80 % power): **4.17 % NO2 per pp**. Estimates below that size can't be told apart from zero by this design; a null here is not evidence of no effect.
Sites: 18, years 2019-2023 (Jan-Apr). Spatial-block bootstrap: 500 draws of 4x4-pixel blocks within site.

| outcome                                            | model                                            | term                         |   estimate |   boot_se |   ci95_lo |   ci95_hi |   n_draws |
|:---------------------------------------------------|:-------------------------------------------------|:-----------------------------|-----------:|----------:|----------:|----------:|----------:|
| raw Jan-Apr NO2 (log of seasonal mean)             | panel two-way FE + SLX                           | L (own built-up)             |     -0.772 |     0.613 |    -1.78  |     0.613 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | panel two-way FE + SLX                           | WL (neighbour built-up)      |      0.835 |     1.38  |    -1.466 |     4.018 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2021-2022 outcomes)                    | L                            |      0.005 |     1.24  |    -2.282 |     2.693 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2021-2022 outcomes)                    | L_lead (placebo: must be ~0) |     -3.564 |     1.6   |    -6.487 |    -0.431 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2021-2022 outcomes)                    | L_lag                        |      0.153 |     1.261 |    -2.815 |     2.013 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2021-2022 outcomes)                    | WL                           |     -3.739 |     3.391 |   -10.462 |     2.575 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | long difference 2019->2023 + site FE (G1 design) | dL                           |      0.131 |     0.634 |    -1.048 |     1.466 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | panel two-way FE + SLX                           | L (own built-up)             |     -0.655 |     0.563 |    -1.533 |     0.608 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | panel two-way FE + SLX                           | WL (neighbour built-up)      |      2.006 |     1.114 |    -0.119 |     4.228 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2021-2022 outcomes)                    | L                            |      0.045 |     1.06  |    -2.147 |     1.99  |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2021-2022 outcomes)                    | L_lead (placebo: must be ~0) |     -4.252 |     1.501 |    -7.529 |    -1.583 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2021-2022 outcomes)                    | L_lag                        |      0.974 |     1.166 |    -1.543 |     3.295 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2021-2022 outcomes)                    | WL                           |     -2.168 |     3.427 |    -8.842 |     4.301 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | long difference 2019->2023 + site FE (G1 design) | dL                           |      0.111 |     0.677 |    -1.17  |     1.415 |       500 |

- raw Jan-Apr NO2 (log of seasonal mean): spatial placebo (500 draws) median b = -0.030, 95 % of |b| below 1.165; permutation p for the panel b = 0.220.
- ERA5-deweathered Jan-Apr NO2 (mean of window logs): spatial placebo (500 draws) median b = 0.057, 95 % of |b| below 1.183; permutation p for the panel b = 0.272.
