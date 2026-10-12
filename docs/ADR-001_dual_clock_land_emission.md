# ADR-001: Dual-clock design for learning land-use change → emission change

**Status:** Proposed
**Date:** 2026-10-12
**Deciders:** Naresh (owner), project guide/examiner
**Supersedes:** "Stage B uses S2 as a static map" in `plans/restart_plan_simple.md`. Step 1 data prep and Step 3 metrics there still apply.

---

## 1. Context

### What Naresh wants
The primary goal is no longer "forecast the next 5-day gas map". It is:

> Learn how **land transformation** (urban expansion, roads, industrial commissioning, crop rotation) drives **emission shifts**, and use that to improve multi-horizon forecasts.

Sentinel-2 must not be a static background map, and land change must not be squeezed into a 10–20 day ConvGRU window.

### What we already know (measured, 2026-10-08)
| Fact | Consequence for this design |
|---|---|
| 12×12 grid at ~5.5 km, 438 windows, ~290 train windows | Very small sample count for any model trained on S5P maps |
| No significant 2019–2024 trend in city-mean NO2/CO/SO2 | Any land-driven signal is **local** (some pixels), not city-wide |
| Damped persistence beats climatology at 5 d (NO2 +6.6%, CO +9.4%); gone by ~30 d / ~5 d | The short-range ceiling without weather is low |
| Learned long-range models failed on 2024 test (NO2 Ridge −4.7%) | Long-range gains need a real new signal, not a bigger model |
| Static S2 explains *where* NO2 sits (NDVI rho −0.53), not *when* | Static land ≠ land **change**. Change was never tested. |
| 2024 test already used once | New clean test = 2025 S5P |

### The dual-clock idea is correct
Land changes over months to years. Gas anomalies change over hours to days. Mixing both into one short window makes the network ignore S2, because cloud and illumination noise in short windows is far larger than real structural change. **This ADR keeps that insight and changes how it is implemented.**

---

## 2. Review of the proposed v2 pipeline

### Keep
1. Two clocks, two streams, separate time resolutions.
2. Multi-temporal, cloud-masked S2 composites, not per-window S2.
3. Physical change indices (ΔNDVI, ΔNDBI, ΔBSI, ΔMNDWI).
4. ERA5 meteorology. It is the biggest missing driver.
5. Climatology + residual formulation.
6. Lead-conditioned blending between short-range and long-range parts.
7. No coordinate channels; test for location memorisation.

### Fix or drop

| # | Proposal item | Problem | Why it matters | Change |
|---|---|---|---|---|
| 1 | Raw R² targets 0.75–0.85 at 5 d | Seasonal cycle alone gives high raw R² (CO is ~83% seasonal). Measured skill vs climatology is +6.6%. | Targets look good and mean nothing; the examiner will ask "vs what?" | Main metric stays **skill vs climatology**. Targets in §6 are hypotheses. |
| 2 | "+0.05–0.10 R² over climatology at 120 d from land change" | No city-scale trend exists. Learned long-range models already failed. | Promising this sets up a negative result as a failure. | Test it as hypothesis H3, stratified by how much a pixel changed. |
| 3 | Macro stream answers the question inside one end-to-end net | Land change → emission is identified from **144 pixels × ~6 years**, heavily spatially correlated. Effective N is tens, not thousands. A ResNet18 + cross-attention trunk cannot learn a causal-looking link from that; it will memorise pixel identity. | This is the core scientific claim. It needs a design where the estimate is defensible. | Answer it with a **panel model (ML/statistics)** on deweathered NO2. Use DL where data is large: S2 at 10 m. |
| 4 | ERA5 as 10 spatial channels on a 12×12 grid | ERA5 is 0.25° (~28 km). The domain is ~3×3 ERA5 cells. The "maps" are nearly uniform. | Wasted channels; the CNN cannot learn spatial advection from them. | Treat ERA5 as **city-scale scalars**: FiLM conditioning for the CNN, plain features for GBM. |
| 5 | Using concurrent ERA5 in a forecast | ERA5 at the target window is not known at issue time. | Silent leakage. Inflates 5–10 d skill. | Two modes, labelled: **(a) deweathering/attribution** may use concurrent ERA5; **(b) forecast** uses only past ERA5, or ECMWF forecasts if obtained. |
| 6 | `L_advect = |u·∂C/∂x + v·∂C/∂y|` | This is not continuity. It pushes gradients to zero or perpendicular to wind. At 4 m/s, air crosses the 66 km domain in ~5 h, so a 5-day composite holds no advection memory. | Wrong physics with confident naming. | **Drop.** If physics is wanted, use the flux-divergence emission estimate (§4, Phase 1 stretch) on daily data. |
| 7 | `L_attrib` = cosine(vision embedding, emission embedding) | No negatives: both encoders collapse to a constant vector and the loss hits 0. | Learns nothing, reports a perfect number. | **Drop.** The panel model in §4.3 answers the attribution question directly. Optional: InfoNCE with negatives, as an ablation only. |
| 8 | NLL with Charbonnier residual | `charbonnier_res ≈ |r|` divided by variance is neither Gaussian nor Laplace NLL. The empty-mask branch returns a tensor, not the tuple. | Miscalibrated uncertainty; crash on all-cloud batch. | Use quantile (pinball) loss for 10/50/90%, which matches the calibrated-band work already done. |
| 9 | Five losses with hand-set lead weights | ~20 tunable weights on ~290 maps tuned on one validation year. | Overfits to 2023, exactly what happened to the Ridge long-range model. | Masked pinball loss only. Add one extra term only if an ablation shows a gain. |
| 10 | Cross-attention fusion on 12×12 | 144 tokens, tiny data. Attention adds parameters without a reason. | Overfit risk, no interpretability gain. | Lead-gated sum (already in the proposal) or FiLM. |
| 11 | "Translation-invariant convs = no memorisation" | Zero padding leaks absolute position on a 12×12 grid. Each pixel's S2 signature is unique, so it acts as an ID. | The claim is untested by design. | Prove it by **spatial-block CV** and, ideally, **leave-one-city-out** (§4.5). |
| 12 | CDS request for 6 years in one call, `format` key | Exceeds new CDS request limits; `format` is replaced by `data_format`/`download_format`. Box has no margin. | Download fails. | Request per year; pad the box by 1° for upwind context. |
| 13 | CO and SO2 in the land→emission link | CO lives 1–2 months and mixes regionally. SO2 is noise-dominated at 5.5 km. | Neither can show a local land signal. | Land→emission study is **NO2 only**. CO/SO2 keep the frozen forecast models. |

### Two hidden confounders the proposal misses (high priority)
These create **fake land-use change or fake emission change**. Both would produce a "significant" result that is an artefact.

1. **TROPOMI NO2 processor versions.** The operational OFFL record mixes versions; v1.4 (Dec 2020) and v2.x changed tropospheric columns, more in polluted and cloudy scenes. Those scenes correlate with built-up land. A step change in 2020–2021 would look like "urban pixels got dirtier".
   **Fix:** use a consistently reprocessed record (S5P-PAL / RPRO v2.4+) for the multi-year work, or at minimum add a version indicator and show the result survives it. First check which versions your GEE exports came from.
2. **Sentinel-2 processing baseline 04.00 (Jan 2022)** added a +1000 offset to L2A reflectance. Raw S2_SR composites show a fake jump in every index in 2022.
   **Fix:** use `COPERNICUS/S2_SR_HARMONIZED`, and build composites from the **same season every year** (e.g. Jan–Mar dry season), so phenology and monsoon greening do not look like land change.

---

## 3. Decision

Split the problem into **three components with different data regimes**, and use each method where it has enough data:

```
                SLOW CLOCK (annual)                              FAST CLOCK (5-day)
 ┌────────────────────────────────────────────┐   ┌───────────────────────────────────────┐
 │ C1  Land-change encoder        [DL]        │   │ C3  Short-range anomaly model [ML→DL] │
 │  S2 10 m, same-season annual composites    │   │  past S5P anomalies + masks           │
 │  → land fractions + change per S5P pixel   │   │  + past ERA5 scalars (+ECMWF fcst)    │
 │    per year  L[p, y]  (millions of pixels) │   │  → A_hat[p, t+h]                      │
 └──────────────────────┬─────────────────────┘   └──────────────────┬────────────────────┘
                        ▼                                            │
 ┌────────────────────────────────────────────┐                      │
 │ C2  Land → emission panel model  [ML/stats]│                      │
 │  deweathered NO2 annual mean per pixel     │                      │
 │  ~ ΔL[p,y] + neighbour ΔL + pixel FE       │                      │
 │    + year FE                               │                      │
 │  → β (effect size) and baseline shift      │                      │
 │    B[p, y]                                 │                      │
 └──────────────────────┬─────────────────────┘                      │
                        ▼                                            ▼
     Forecast(p, t+h) = Clim[p, season(t+h)] + B[p, year(t+h)] + g(h) · A_hat[p, t+h]
                         static seasonal      structural shift    fast anomaly, fades with lead
```

Why this split:
- **DL where data is big:** S2 at 10 m gives millions of labelled pixels per year (labels from Dynamic World / Open Buildings Temporal). Segmentation and change detection are what CNNs and foundation models are good at.
- **ML/statistics where data is small:** the land→emission link has tens of effective samples. A two-way fixed-effects panel model removes the static "where" confound (pixel FE) and city-wide shocks like COVID 2020, economy and retrieval steps (year FE). What is left is **within-pixel change vs within-pixel change**. That is the defensible version of "land change drives emission change".
- **Deweathering first:** weather moves NO2 by far more than one year of land change does. Removing weather-driven variance with ERA5 + GBM is the step that makes the slow signal visible. It is a standard, citable method (meteorological normalisation, Grange et al. 2018; `rmweather`).
- **Forecasting uses the result:** the structural shift `B` replaces the frozen "climatology never changes" assumption. That is exactly where land change can help long leads, and only there.

---

## 4. Components

### 4.1 Data (Phase 1)

| Data | Source | Clock | Notes |
|---|---|---|---|
| S5P NO2 tropospheric column | Reprocessed record (PAL/RPRO) if possible; else OFFL with version flag | 5-day (fast), annual (slow) | Extend to Jul 2018 → 2025. 2025 = clean test. |
| S5P daily L3 NO2 (stretch) | GEE `COPERNICUS/S5P/OFFL/L3_NO2` | daily | Only for flux-divergence emissions (§4.2) |
| S2 L2A | `COPERNICUS/S2_SR_HARMONIZED`, SCL cloud mask | annual, same season | 2018/19 → 2025 |
| Annual S2 embeddings (strong option) | `GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL` (AlphaEarth, 64-d, 10 m) | annual | Year-comparable by design. Change = embedding difference. |
| Building presence/height per year | `GOOGLE/Research/open-buildings-temporal/v1` | annual (2016–2023) | Direct urban-growth label; check coverage of your box |
| Land-cover probabilities | Dynamic World (per scene → annual mean probability) | annual | Labels for C1 |
| Night lights | VIIRS DNB annual/monthly | annual | Activity proxy; complements reflectance |
| Meteorology | ERA5 single levels: u10, v10, BLH, t2m, tp, sp; plus 850 hPa u/v | hourly → overpass (~08 UTC) and 5-day means | Box padded by 1°; one CDS request per year |

ERA5 features (city scalars per window): wind speed, wind direction as sin/cos, BLH at overpass, t2m, precipitation sum, ventilation index = BLH × wind speed, count of rainy days. These go in as scalars, not maps.

### 4.2 Deweathering (Phase 2, ML)
Goal: per pixel, a weather-normalised NO2 series where remaining change reflects emissions.

- Model per pixel (or pooled with pixel ID as a categorical): LightGBM `NO2[p,t] ~ ERA5[t] + day-of-year + time index`.
- Normalise: re-predict each window many times with ERA5 resampled from the same season across all years; average. Output `NO2_dw[p, t]`.
- Aggregate to annual (and dry-season) means `Y[p, y]`.
- Check: deweathered series should show the 2020 lockdown dip **more** clearly than raw. If it does not, the deweathering is not working.

Stretch (only if Phases 2–3 show a signal): **flux-divergence emission estimate** (Beirle et al. 2019/2021) from daily L3 NO2 and ERA5 winds, `E ≈ ∇·(V·u) + V/τ`. That converts columns (smeared by transport) into emissions, which is what land change actually drives. This is the correct place for wind physics, not a CNN loss term.

### 4.3 C1: land-change encoder (Phase 3, DL)
Two tiers. Do Tier 1 first; Tier 2 only if time allows.

**Tier 1 (cheap, strong):** frozen annual AlphaEarth embeddings, plus Dynamic World annual probabilities, plus Open Buildings Temporal built-up fraction. Aggregate each 10 m layer into each S5P pixel (area mean) per year, giving `L[p, y]`, and changes `ΔL[p, y] = L[p, y] − L[p, 2019]`.

**Tier 2 (your own DL contribution):** a small Siamese U-Net on bi-temporal same-season S2 composites (bands B2, B3, B4, B8, B11, B12 + indices), trained to predict built-up / vegetation / bare change, with Dynamic World or Open Buildings Temporal change as labels.
- Train on 64×64 or 128×128 patches; there are millions of 10 m pixels, so DL is justified here.
- Initialise from an EO foundation model (SSL4EO-S12 ResNet or Prithvi-EO-2.0) if it helps; compare against training from scratch.
- Score with a **held-out spatial block**: per-class F1, change-F1, kappa.
- Output per S5P pixel per year: % built-up, % vegetation, % bare, % new built-up since reference year.

C1 is a strong result on its own, independent of whether C2 finds an effect.

### 4.4 C2: land → emission panel model (Phase 4, the primary scientific result)
```
Y[p,y] = α_p + γ_y + β · ΔL[p,y] + θ · W·ΔL[p,y] + δ · NTL[p,y] + ε[p,y]
```
- `α_p` pixel fixed effect: removes static "where" (greenness, coast, old industry).
- `γ_y` year fixed effect: removes COVID 2020, economy, city-wide meteorology left after deweathering, and the city-wide part of retrieval-version steps.
- `W·ΔL` = land change in the ring of neighbouring pixels. NO2 at a pixel is partly emitted upwind; this is a spatial-lag (spatial Durbin) term.
- Standard errors: cluster by spatial block, or use spatial-block bootstrap. Naive SEs will be far too small.
- Non-linear check: LightGBM on the same within-transformed variables plus SHAP. It must agree in sign with the linear model.

Validation that the effect is real:
1. **Placebo:** replace ΔL with land change shuffled between pixels. β should vanish.
2. **Lead/lag test:** future land change must not predict past NO2. If it does, the link is confounded.
3. **Event studies:** 2–4 sites with documented change inside the box (e.g. new industrial/logistics parks in the Sriperumbudur–Oragadam belt, North Chennai/Ennore power and port expansion). Plot deweathered NO2 before and after the dated change against matched stable pixels. Verify commissioning dates from public sources before using them.
4. **Power check (do this first, in Gate G1):** if most pixels have near-zero ΔL, there is nothing to learn. Report the distribution of ΔL before fitting anything.

### 4.5 Fixing the sample-size problem: multi-city panel (gated)
If Gate G1 shows weak power in Chennai alone, add 4–6 Indian cities with similar sensors and very different growth (e.g. Bengaluru, Hyderabad, Pune, Ahmedabad, Kolkata). The pipeline is identical; data at 5.5 km is small.
- Effective N grows by roughly the number of cities.
- **Leave-one-city-out** evaluation is the real proof of "no coordinate memorisation" and "learns the archetype of an emitter, not the location".
- Chennai stays the case study and the forecast target.

### 4.6 C3: short-range anomaly model (Phase 5)
Ladder unchanged; climb only if the level below is beaten on 2023 validation **and** 2025 test.

| Level | Model | Inputs |
|---|---|---|
| 0 | Climatology | — |
| 1 | Damped persistence (frozen champion) | `A[t]` |
| 2 | LightGBM / Ridge | lags, 3×3 neighbour mean, city mean, masks, **past ERA5 scalars**, `L[p,y]` |
| 3 | Small U-Net/ConvGRU with FiLM | anomaly maps + masks as channels; ERA5 scalars and lead `h` via FiLM; `L[p,y]` maps as static channels |

For Level 3, the DL advantage that is actually available: train on **daily** L3 maps (≈5× more samples than 5-day composites), with masks, to predict the next 5-day composite. Without that, expect it to tie Level 2.

If ECMWF medium-range forecasts can be obtained for 2019–2025, add them as a labelled "forecast-met" variant. Report it separately from the past-only variant.

### 4.7 Fusion and multi-horizon output (Phase 6)
```
F[p, t+h] = Clim[p, s(t+h)] + B[p, y(t+h)] + g(h) · A_hat[p, t+h]
```
- `B[p,y] = β̂ · ΔL[p,y] + θ̂ · W·ΔL[p,y]` from C2, using the latest observed land state at issue time (no future land).
- `g(h)` is learned per lead on validation (monotone, 0 ≤ g ≤ 1). This replaces the proposal's hand-set gate.
- Bands: quantile models per lead, as already done; target 80% coverage.
- Horizons: 5, 10, 20, 30, 60, 120, 150 d (h = 1…30 windows).

---

## 5. Options considered

### Option A: proposed v2 end-to-end dual-branch network
| Dimension | Assessment |
|---|---|
| Complexity | High (two encoders, cross-attention, 5 losses, ~20 weights) |
| Data fit | Poor: ~290 maps, effective N tens for the land link |
| Scientific defensibility | Low: no control for static confound, year shocks, retrieval version |
| Interpretability | Low: effect of land change is buried in a latent space |

**Pros:** one model, looks advanced, matches the dual-clock story.
**Cons:** cannot separate land effect from location memorisation; two broken loss terms; leakage risk with ERA5; promised skill unsupported by measured results.

### Option B (chosen): three components, DL on S2, panel ML for the link, fusion for forecasting
| Dimension | Assessment |
|---|---|
| Complexity | Medium; each part testable alone |
| Data fit | Good: each method sits where it has enough data |
| Scientific defensibility | High: fixed effects, placebo, lead/lag, event studies, block CV |
| Interpretability | High: β is an effect size you can state in one sentence |

**Pros:** primary goal answered directly; negative result is still publishable; DL used where it wins.
**Cons:** less "one big network"; needs deweathering and reprocessed S5P; land effect may be small.

### Option C: keep S2 static, add ERA5 only
| Dimension | Assessment |
|---|---|
| Complexity | Low |
| Data fit | Good |
| Scientific defensibility | Good for short range |
| Goal fit | **Fails the primary goal** (no land change) |

Kept as the ablation baseline inside Option B.

---

## 6. Hypotheses and honest targets

All on skill vs climatology unless stated. Ranges are guesses to test, not promises.

| ID | Hypothesis | Test | Plausible outcome |
|---|---|---|---|
| H1 | Deweathering exposes local NO2 change hidden by weather | Variance of annual pixel means drops; lockdown dip sharper | Likely |
| H2 | Within-pixel land change predicts within-pixel deweathered NO2 change (β > 0 for new built-up / industrial) | Panel FE + placebo + lead/lag + events | Uncertain; strongest in fast-changing peri-urban pixels; may be null in Chennai alone |
| H3 | Structural shift `B` improves long-lead forecasts | 60–150 d skill vs static climatology, **split by pixel change tercile** | ~0 city-wide; small positive only in high-change pixels, if H2 holds |
| H4 | Past ERA5 improves 5 d NO2 skill | Level 2 with vs without ERA5 | Modest gain over +6.6% |
| H5 | Forecast meteorology helps 5–10 d | Forecast-met variant | Larger gain, if data obtainable |

A null H2 with proper controls is a valid, defensible thesis result: "at 5.5 km and six years, land change effects on NO2 are below detection after controlling for weather and year shocks; here is the detection limit".

---

## 7. Evaluation protocol
- Splits: train 2018-07–2022, validate 2023–2024, **test 2025** (2024 was seen once).
- C1: spatial-block holdout; F1, change-F1, kappa.
- C2: spatial-block bootstrap CIs; placebo; lead/lag; leave-one-city-out if multi-city.
- C3 and fusion: skill vs climatology per lead; CRPS; band coverage; per-season (monsoon vs dry); **stratified by land-change tercile**.
- Ablations: S2 static vs S2 change; with/without ERA5; with/without `B`; retrieval-version flag on/off; exclude 2020 Q2.
- Freeze, retrain on train+val, score 2025 once.

---

## 8. Consequences
- **Easier:** answering "does land change drive emissions?" with one defensible number; explaining results to an examiner; reusing the frozen forecast models as baselines.
- **Harder:** data work (reprocessed S5P, harmonised S2, ERA5, deweathering) before any model.
- **Revisit:** multi-city extension after Gate G1; flux-divergence emissions after Gate G2; daily-data CNN after Level 2 results.

---

## 9. Action items (gated)

| Phase | Tasks | Exit gate |
|---|---|---|
| 1 Data | Check S5P versions in existing exports; get PAL/RPRO NO2 or add version flag. Same-season harmonised S2 composites 2019–2025. ERA5 per year (padded box). AlphaEarth annual embeddings, Dynamic World annual probabilities, Open Buildings Temporal, VIIRS. Download 2025 S5P. | All layers on the 12×12 grid; versions documented |
| **G1 Power** | Distribution of ΔL per pixel; how many pixels changed > X% built-up. | If < ~10 pixels changed meaningfully → add cities (§4.5) before C2 |
| 2 Deweather | LightGBM meteorological normalisation per pixel; annual `Y[p,y]`. | Lockdown dip sharper than raw |
| 3 C1 | Tier 1 features; Tier 2 Siamese U-Net if time. | Spatial-block F1 reported |
| 4 C2 | Two-way FE panel + spatial lag; placebo; lead/lag; 2–4 event studies. | β with block-bootstrap CI; placebo null |
| **G2** | Is β significant and sign-consistent? | Yes → use `B` in fusion. No → `B = 0`, report detection limit |
| 5 C3 | Level 2 with past ERA5; Level 3 only if Level 2 beats damped persistence. | Beats Level 1 on val and test |
| 6 Fusion | Learn `g(h)`; quantile bands; skill-vs-lead by change tercile. | Coverage ~80% all leads |
| 7 Report | Single 2025 test; figures; write-up. | — |

Recommended order if time is short: **Phase 1 → G1 → Phase 2 → Phase 4 → Phase 3 Tier 1 → Phase 5 → Phase 6.** Tier 2 DL change detection is the first thing to cut, and the first thing to add back if there is time.

---

## References to cite
- Grange, Carslaw et al. (2018), random-forest meteorological normalisation; `rmweather`.
- Beirle et al. (2019, Sci. Adv.; 2021, ESSD), flux-divergence NOx emissions from TROPOMI.
- van Geffen et al. (2022, AMT), TROPOMI NO2 v2.2 changes; S5P-PAL reprocessing docs.
- Brown et al. (2022, Sci. Data), Dynamic World.
- Google Open Buildings 2.5D Temporal; AlphaEarth Foundations / Satellite Embedding V1.
- Islam et al. (2020, ICLR), position information encoded by CNN padding.
