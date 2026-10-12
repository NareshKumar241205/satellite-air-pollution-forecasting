# C2 panel: EXPLORATORY (G1 failed at 0.5 %/pp; MDE 0.73 %/pp); raw and ERA5-deweathered Jan-Apr NO2

G1 detection limit (pooled long difference, 80 % power): **4.17 % NO2 per pp**. Estimates below that size can't be told apart from zero by this design; a null here is not evidence of no effect.
Sites: 18, years 2019-2023 (Jan-Apr). Spatial-block bootstrap: 500 draws of 4x4-pixel blocks within site.

| outcome                                            | model                                            | term                         |   estimate |   boot_se |   ci95_lo |   ci95_hi |   n_draws |
|:---------------------------------------------------|:-------------------------------------------------|:-----------------------------|-----------:|----------:|----------:|----------:|----------:|
| raw Jan-Apr NO2 (log of seasonal mean)             | panel two-way FE + SLX                           | L (own built-up)             |     -0.388 |     0.538 |    -1.259 |     0.871 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | panel two-way FE + SLX                           | WL (neighbour built-up)      |      1.662 |     1.294 |    -0.561 |     4.4   |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2020-2022 outcomes)                    | L                            |     -0.229 |     0.661 |    -1.552 |     1.102 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2020-2022 outcomes)                    | L_lead (placebo: must be ~0) |      2.745 |     1.208 |     0.357 |     4.868 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2020-2022 outcomes)                    | L_lag                        |      5.089 |     0.933 |     3.373 |     6.835 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2020-2022 outcomes)                    | WL                           |      3.224 |     2.429 |    -1.259 |     8.118 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | long difference 2019->2023 + site FE (G1 design) | dL                           |      0.131 |     0.634 |    -1.048 |     1.466 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | panel two-way FE + SLX                           | L (own built-up)             |     -0.283 |     0.522 |    -1.137 |     0.883 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | panel two-way FE + SLX                           | WL (neighbour built-up)      |      3.025 |     1.115 |     0.885 |     5.393 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2020-2022 outcomes)                    | L                            |     -0.18  |     0.452 |    -1.083 |     0.741 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2020-2022 outcomes)                    | L_lead (placebo: must be ~0) |      0.602 |     1.294 |    -2.36  |     2.701 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2020-2022 outcomes)                    | L_lag                        |      6.629 |     1.15  |     4.549 |     8.983 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2020-2022 outcomes)                    | WL                           |      7.539 |     2.437 |     2.889 |    12.456 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | long difference 2019->2023 + site FE (G1 design) | dL                           |      0.111 |     0.677 |    -1.17  |     1.415 |       500 |

- raw Jan-Apr NO2 (log of seasonal mean): spatial placebo (500 draws) median b = 0.012, 95 % of |b| below 1.058; permutation p for the panel b = 0.502.
- ERA5-deweathered Jan-Apr NO2 (mean of window logs): spatial placebo (500 draws) median b = 0.025, 95 % of |b| below 0.981; permutation p for the panel b = 0.550.
