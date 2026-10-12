# C2 panel: EXPLORATORY (G1 failed at 0.5 %/pp; MDE 0.73 %/pp); raw and ERA5-deweathered Jan-Apr NO2

**Units:** the land variable is 10·ln(1 + night-light radiance), so b is % NO2 per 10 % change in lights, not per percentage point of built-up area. The detection limit is not estimable under the G1 noise rule (see FINAL_REPORT §4c); an optimistic bound from the unadjusted spread is 0.77 % per 10 % light.

G1 detection limit (pooled long difference, 80 % power): **not estimable (no real land signal after noise removal) % NO2 per pp**. Estimates below that size can't be told apart from zero by this design; a null here is not evidence of no effect.
Sites: 18, years 2019-2023 (Jan-Apr). Spatial-block bootstrap: 500 draws of 4x4-pixel blocks within site.

| outcome                                            | model                                            | term                         |   estimate |   boot_se |   ci95_lo |   ci95_hi |   n_draws |
|:---------------------------------------------------|:-------------------------------------------------|:-----------------------------|-----------:|----------:|----------:|----------:|----------:|
| raw Jan-Apr NO2 (log of seasonal mean)             | panel two-way FE + SLX                           | L (own built-up)             |      0.106 |     0.074 |    -0.037 |     0.25  |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | panel two-way FE + SLX                           | WL (neighbour built-up)      |      0.926 |     0.314 |     0.291 |     1.52  |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2020-2022 outcomes)                    | L                            |      0.242 |     0.155 |    -0.06  |     0.539 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2020-2022 outcomes)                    | L_lead (placebo: must be ~0) |      0.253 |     0.217 |    -0.249 |     0.593 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2020-2022 outcomes)                    | L_lag                        |     -0.526 |     0.227 |    -0.867 |    -0.017 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2020-2022 outcomes)                    | WL                           |      1.131 |     0.416 |     0.376 |     1.906 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | long difference 2019->2023 + site FE (G1 design) | dL                           |      0.079 |     0.079 |    -0.076 |     0.232 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | panel two-way FE + SLX                           | L (own built-up)             |      0.085 |     0.079 |    -0.065 |     0.242 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | panel two-way FE + SLX                           | WL (neighbour built-up)      |      0.823 |     0.272 |     0.278 |     1.323 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2020-2022 outcomes)                    | L                            |      0.107 |     0.137 |    -0.161 |     0.385 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2020-2022 outcomes)                    | L_lead (placebo: must be ~0) |      0.116 |     0.206 |    -0.321 |     0.481 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2020-2022 outcomes)                    | L_lag                        |     -0.184 |     0.211 |    -0.552 |     0.233 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2020-2022 outcomes)                    | WL                           |      0.884 |     0.36  |     0.188 |     1.57  |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | long difference 2019->2023 + site FE (G1 design) | dL                           |      0.084 |     0.093 |    -0.095 |     0.267 |       500 |

- raw Jan-Apr NO2 (log of seasonal mean): spatial placebo (500 draws) median b = 0.021, 95 % of |b| below 0.152; permutation p for the panel b = 0.160.
- ERA5-deweathered Jan-Apr NO2 (mean of window logs): spatial placebo (500 draws) median b = 0.013, 95 % of |b| below 0.149; permutation p for the panel b = 0.312.
