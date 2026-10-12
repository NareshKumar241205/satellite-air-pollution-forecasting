# Lead-test diagnostics (plan: `Land Change/DIAGNOSTICS_PLAN.md`)

## D1: scene/atmosphere covariates (within site-year r)

| covariate       |      L |   L_lead |   dL_next |   y_res |
|:----------------|-------:|---------:|----------:|--------:|
| aod055          |  0.556 |    0.593 |     0.096 |   0.018 |
| aod055_next     |  0.534 |    0.554 |     0.103 |  -0.02  |
| brightness      |  0.506 |    0.498 |     0.102 |   0.021 |
| brightness_next |  0.53  |    0.536 |     0.102 |  -0.004 |
| cloud_prob      | -0.022 |   -0.03  |     0.006 |   0.017 |
| cloud_prob_next | -0.014 |   -0.024 |    -0.033 |   0.007 |
| dw_images       | -0.086 |   -0.088 |    -0.029 |  -0.008 |
| dw_images_next  | -0.089 |   -0.089 |    -0.014 |  -0.007 |

## D2: level dependence of the yearly label change (within site)

| pair       |   dL_per_100logNO2 |   dL_on_L0 |   dL_sd_within |    n |
|:-----------|-------------------:|-----------:|---------------:|-----:|
| 2019->2020 |              0.804 |      0.005 |          1.429 | 2227 |
| 2020->2021 |              1.219 |      0.003 |          1.388 | 2227 |
| 2021->2022 |              1.975 |      0.032 |          1.739 | 2227 |
| 2022->2023 |             -0.658 |     -0.027 |          2.023 | 2227 |

## D3: pixel-specific linear trends

`{"lead_with_pixel_trends": 0.11, "L": -0.015, "n": 6681, "note": "outcomes 2020-2022: 3 periods per pixel, so 1 residual df per pixel after a pixel trend"}`

## D4: leave one out (lead coefficient)

| dropped              |   lead |   change_vs_full_pct |
|:---------------------|-------:|---------------------:|
| none                 |  0.308 |                0     |
| site Chennai         |  0.381 |               23.566 |
| site Coimbatore      |  0.321 |                4.282 |
| site Cuddalore       |  0.332 |                7.583 |
| site Dindigul        |  0.312 |                1.308 |
| site Hosur           |  0.3   |               -2.578 |
| site Kanchipuram     |  0.277 |              -10.225 |
| site Karur           |  0.316 |                2.657 |
| site Madurai         |  0.297 |               -3.729 |
| site Nagapattinam    |  0.295 |               -4.335 |
| site Ramanathapuram  |  0.326 |                5.741 |
| site Salem           |  0.236 |              -23.311 |
| site Sivakasi        |  0.321 |                4.302 |
| site Thanjavur       |  0.345 |               11.781 |
| site Thoothukudi     |  0.306 |               -0.583 |
| site Tiruchirappalli |  0.3   |               -2.525 |
| site Tiruppur        |  0.293 |               -5.009 |
| site Tiruvannamalai  |  0.3   |               -2.566 |
| site Vellore         |  0.281 |               -8.771 |
| outcome year 2020    |  0.009 |              -97.1   |
| outcome year 2021    |  0.343 |               11.357 |
| outcome year 2022    |  0.464 |               50.636 |

## D5: drop the outer ring of blocks

`{"lead_full": 0.308, "lead_inner": 0.273, "change_pct": -11.373}`

## D6: residual spatial structure

`{"e_y~e_y": 0.724, "e_lead~e_lead": 0.483, "e_y~e_lead(neighbour)": 0.023, "e_y~e_lead(same pixel)": 0.0}`
