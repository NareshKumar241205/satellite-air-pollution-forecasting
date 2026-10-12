# C2 panel: POST-HOC DIAGNOSTIC D7, outcome years [2020] removed (DIAGNOSTICS_PLAN.md addendum). EXPLORATORY (G1 failed at 0.5 %/pp; MDE 0.73 %/pp); raw and ERA5-deweathered Jan-Apr NO2

G1 detection limit (pooled long difference, 80 % power): **0.73 % NO2 per pp**. Estimates below that size can't be told apart from zero by this design; a null here is not evidence of no effect.
Sites: 18, years 2019-2023 (Jan-Apr). Spatial-block bootstrap: 500 draws of 4x4-pixel blocks within site.

| outcome                                            | model                                            | term                         |   estimate |   boot_se |   ci95_lo |   ci95_hi |   n_draws |
|:---------------------------------------------------|:-------------------------------------------------|:-----------------------------|-----------:|----------:|----------:|----------:|----------:|
| raw Jan-Apr NO2 (log of seasonal mean)             | panel two-way FE + SLX                           | L (own built-up)             |      0.001 |     0.044 |    -0.078 |     0.099 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | panel two-way FE + SLX                           | WL (neighbour built-up)      |      0.254 |     0.118 |     0.024 |     0.487 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2021-2022 outcomes)                    | L                            |     -0.029 |     0.108 |    -0.24  |     0.171 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2021-2022 outcomes)                    | L_lead (placebo: must be ~0) |      0.009 |     0.132 |    -0.27  |     0.247 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2021-2022 outcomes)                    | L_lag                        |     -0.246 |     0.157 |    -0.558 |     0.055 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | lead/lag (2021-2022 outcomes)                    | WL                           |     -0.107 |     0.195 |    -0.527 |     0.213 |       500 |
| raw Jan-Apr NO2 (log of seasonal mean)             | long difference 2019->2023 + site FE (G1 design) | dL                           |     -0.026 |     0.064 |    -0.151 |     0.102 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | panel two-way FE + SLX                           | L (own built-up)             |      0.011 |     0.044 |    -0.073 |     0.09  |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | panel two-way FE + SLX                           | WL (neighbour built-up)      |      0.351 |     0.104 |     0.142 |     0.544 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2021-2022 outcomes)                    | L                            |     -0.091 |     0.095 |    -0.279 |     0.095 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2021-2022 outcomes)                    | L_lead (placebo: must be ~0) |      0.01  |     0.112 |    -0.21  |     0.239 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2021-2022 outcomes)                    | L_lag                        |     -0.249 |     0.144 |    -0.507 |     0.042 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | lead/lag (2021-2022 outcomes)                    | WL                           |     -0.184 |     0.177 |    -0.519 |     0.145 |       500 |
| ERA5-deweathered Jan-Apr NO2 (mean of window logs) | long difference 2019->2023 + site FE (G1 design) | dL                           |      0.003 |     0.068 |    -0.125 |     0.142 |       500 |

- raw Jan-Apr NO2 (log of seasonal mean): spatial placebo (500 draws) median b = 0.011, 95 % of |b| below 0.079; permutation p for the panel b = 0.984.
- ERA5-deweathered Jan-Apr NO2 (mean of window logs): spatial placebo (500 draws) median b = 0.012, 95 % of |b| below 0.082; permutation p for the panel b = 0.808.
