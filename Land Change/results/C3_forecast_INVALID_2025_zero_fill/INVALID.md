# INVALID: the first 2025 scoring (12 Oct 2026, 22:37)

It was made with a defect in the 2025 NO2 files downloaded by `scripts/ee_s5p_no2_windows.py`:
- Pixels outside the Tamil Nadu clip were filled with **0** instead of NaN.
- Up to ~50 % of the coastal and border boxes (Nagapattinam, Cuddalore, Ramanathapuram, Thoothukudi, Hosur, Chennai) were averaged in as zero NO2.
- That gave apparent 2025 drops of 40–60 %.

The method check missed it, because it compared only the −9999 mask and treated NaN and 0 alike.

The fix:
- the files were repaired into `data/tnaq_extra_v2`;
- the downloader and the check were corrected (the check now compares NaN, −9999 and valid classes, and passes on 2024-01-05).

Nothing in the model or its settings was changed. The 2025 score was re-run once with the repaired data, in
`../C3_forecast/`. These files are kept for the record.
