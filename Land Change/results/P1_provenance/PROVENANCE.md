# TROPOMI NO2 processor versions in the TNAQ export

- Files checked: 473 (2018-07-10 .. 2024-12-30).
- Collection: `COPERNICUS/S5P/OFFL/L3_NO2`, band `tropospheric_NO2_column_number_density` (from `D:\CLG\SEM-7\FINAL PROJECT\AirQ\Scripts\Data_acquisition.js`).
- Distinct metadata sets in the GeoTIFFs: 1. Dataset tags: `[('AREA_OR_POINT', 'Area')]`; bands: `('value', 'n_obs')`.
- Manifest columns: `system:index, end_date, n_images, product, start_date, valid_fraction, .geo`.

**Result: the processor version is not recorded anywhere in the local export.** Neither the GeoTIFFs nor
the manifests carry it. The exact per-window versions need `scripts/ee_s5p_versions.py` (Earth Engine).

The collection is the operational OFFL stream, not a consistent reprocessing, so by the public timeline
(approximate, in `provenance.py`) the Jan-Apr season windows fall into:

|   year |   1.2.x |   1.3.x |   1.4.0 |   2.3.1 |   2.4.0 |   2.5.0 and later (approximate) |
|-------:|--------:|--------:|--------:|--------:|--------:|--------------------------------:|
|   2019 |      16 |       8 |       0 |       0 |       0 |                               0 |
|   2020 |       0 |      24 |       0 |       0 |       0 |                               0 |
|   2021 |       0 |       0 |      24 |       0 |       0 |                               0 |
|   2022 |       0 |       0 |       0 |      24 |       0 |                               0 |
|   2023 |       0 |       0 |       0 |       0 |      14 |                              10 |
|   2024 |       0 |       0 |       0 |       0 |       0 |                              24 |

Every 2019-2021 season is v1.x and every 2022-2024 season is v2.3.1 or later, so the **2021 -> 2022 pair
straddles the v1 -> v2 change** (v2.2, 1 Jul 2021: new cloud and surface-albedo inputs that raise
tropospheric NO2, more over polluted pixels). City-year fixed effects absorb a uniform step but not a
level-dependent one, so G1 is also reported without that pair, and the main analysis must either use a
consistently reprocessed product (S5P-PAL / RPRO v2.4) or model the step as level-dependent.

_Exact versions: not yet queried (`scripts/ee_s5p_versions.py` needs an Earth Engine account)._

Mean same-season dNO2 (%, log change) per year pair and city, all available years. A retrieval
step would show up in every city in the same pair (2020 and 2021 also carry the lockdown and
its rebound):

|              |   Chennai |   Coimbatore |   Hosur |   Madurai |   Salem |   Thoothukudi |   Tiruchirappalli |   Tiruppur |   Vellore |   mean_all_cities |
|:-------------|----------:|-------------:|--------:|----------:|--------:|--------------:|------------------:|-----------:|----------:|------------------:|
| (2019, 2020) |     -22.2 |        -24.5 |   -27.7 |     -18.1 |   -23.4 |         -22.4 |             -13.4 |      -19.8 |     -26.1 |             -21.9 |
| (2020, 2021) |      28.6 |         16.3 |    18.1 |      10.4 |    18.9 |          10.3 |              16   |       15   |      11.4 |              16.1 |
| (2021, 2022) |      -2.6 |         -3.4 |    -4.9 |       4.1 |    -7.4 |           4.4 |              -2.7 |       -1.9 |       0.2 |              -1.6 |
| (2022, 2023) |      11.6 |          7.4 |     7.3 |       1.9 |     6   |          21.9 |              14.5 |       10   |      11.8 |              10.3 |
| (2023, 2024) |     -30.6 |         -9.7 |    -1.9 |     -13.1 |    -4.4 |           1.8 |             -17.4 |      -11.3 |     -15.8 |             -11.4 |
