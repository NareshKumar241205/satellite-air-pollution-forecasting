# C2 panel: EXPLORATORY (G1 failed at 0.5 %/pp; MDE 0.73 %/pp); raw and ERA5-deweathered Jan-Apr NO2

G1 detection limit (pooled long difference, 80 % power): **4.73 % SO2 per pp**. Estimates below that size can't be told apart from zero by this design; a null here is not evidence of no effect.
Sites: 18, years 2019-2023 (Jan-Apr). Spatial-block bootstrap: 500 draws of 4x4-pixel blocks within site.

| outcome                                            | model                                            | term                         |   estimate |   boot_se |   ci95_lo |   ci95_hi |   n_draws |
|:---------------------------------------------------|:-------------------------------------------------|:-----------------------------|-----------:|----------:|----------:|----------:|----------:|
| raw Jan-Apr NO2 (log of seasonal mean)             | panel two-way FE + SLX                           | L (own built-up)             |      1.07  |     0.725 |    -0.362 |     2.511 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | panel two-way FE + SLX                           | WL (neighbour built-up)      |     -1.15  |     1.069 |    -3.54  |     0.813 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2020-2022 outcomes)                    | L                            |      0.168 |     0.997 |    -1.773 |     2.248 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2020-2022 outcomes)                    | L_lead (placebo: must be ~0) |     -0.344 |     1.112 |    -2.661 |     1.46  |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2020-2022 outcomes)                    | L_lag                        |      1.266 |     1.144 |    -0.877 |     3.492 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2020-2022 outcomes)                    | WL                           |      0.732 |     2.166 |    -3.545 |     4.58  |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | long difference 2019->2023 + site FE (G1 design) | dL                           |      1.142 |     0.948 |    -0.773 |     2.932 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | panel two-way FE + SLX                           | L (own built-up)             |      0.257 |     0.301 |    -0.324 |     0.83  |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | panel two-way FE + SLX                           | WL (neighbour built-up)      |     -0.011 |     0.439 |    -0.837 |     0.838 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2020-2022 outcomes)                    | L                            |      0.229 |     0.445 |    -0.662 |     1.067 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2020-2022 outcomes)                    | L_lead (placebo: must be ~0) |      0.133 |     0.415 |    -0.653 |     0.901 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2020-2022 outcomes)                    | L_lag                        |      0.442 |     0.453 |    -0.364 |     1.362 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2020-2022 outcomes)                    | WL                           |     -0.196 |     0.701 |    -1.574 |     1.238 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | long difference 2019->2023 + site FE (G1 design) | dL                           |     -0.142 |     0.53  |    -1.206 |     0.854 |       500 |

- raw Jan-Apr NO2 (log of seasonal mean): spatial placebo (500 draws) median b = 0.207, 95 % of |b| below 1.338; permutation p for the panel b = 0.114.
- ERA5-deweathered Jan-Apr NO2 (mean of window logs): spatial placebo (500 draws) median b = 0.061, 95 % of |b| below 0.617; permutation p for the panel b = 0.354.
