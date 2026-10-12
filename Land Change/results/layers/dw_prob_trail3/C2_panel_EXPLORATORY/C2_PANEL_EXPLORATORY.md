# C2 panel: EXPLORATORY (G1 failed at 0.5 %/pp; MDE 0.73 %/pp); raw and ERA5-deweathered Jan-Apr NO2

G1 detection limit (pooled long difference, 80 % power): **1.52 % NO2 per pp**. Estimates below that size can't be told apart from zero by this design; a null here is not evidence of no effect.
Sites: 18, years 2019-2023 (Jan-Apr). Spatial-block bootstrap: 500 draws of 4x4-pixel blocks within site.

| outcome                                            | model                                            | term                         |   estimate |   boot_se |   ci95_lo |   ci95_hi |   n_draws |
|:---------------------------------------------------|:-------------------------------------------------|:-----------------------------|-----------:|----------:|----------:|----------:|----------:|
| raw Jan-Apr NO2 (log of seasonal mean)             | panel two-way FE + SLX                           | L (own built-up)             |     -0.034 |     0.071 |    -0.157 |     0.107 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | panel two-way FE + SLX                           | WL (neighbour built-up)      |      0.168 |     0.21  |    -0.218 |     0.577 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2020-2022 outcomes)                    | L                            |     -0.31  |     0.155 |    -0.627 |     0.008 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2020-2022 outcomes)                    | L_lead (placebo: must be ~0) |      0.238 |     0.268 |    -0.271 |     0.771 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2020-2022 outcomes)                    | L_lag                        |      0.581 |     0.242 |     0.088 |     1.018 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2020-2022 outcomes)                    | WL                           |     -0.177 |     0.338 |    -0.838 |     0.528 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | long difference 2019->2023 + site FE (G1 design) | dL                           |     -0.065 |     0.101 |    -0.275 |     0.126 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | panel two-way FE + SLX                           | L (own built-up)             |     -0.019 |     0.067 |    -0.155 |     0.106 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | panel two-way FE + SLX                           | WL (neighbour built-up)      |      0.689 |     0.205 |     0.269 |     1.065 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2020-2022 outcomes)                    | L                            |     -0.451 |     0.148 |    -0.74  |    -0.152 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2020-2022 outcomes)                    | L_lead (placebo: must be ~0) |      0.543 |     0.26  |     0.023 |     1.024 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2020-2022 outcomes)                    | L_lag                        |      0.638 |     0.224 |     0.226 |     1.055 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2020-2022 outcomes)                    | WL                           |      0.783 |     0.317 |     0.153 |     1.419 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | long difference 2019->2023 + site FE (G1 design) | dL                           |      0     |     0.107 |    -0.217 |     0.192 |       500 |

- raw Jan-Apr NO2 (log of seasonal mean): spatial placebo (500 draws) median b = -0.011, 95 % of |b| below 0.145; permutation p for the panel b = 0.686.
- ERA5-deweathered Jan-Apr NO2 (mean of window logs): spatial placebo (500 draws) median b = -0.000, 95 % of |b| below 0.143; permutation p for the panel b = 0.786.
