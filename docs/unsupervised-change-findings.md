# Kidepo landscape LCLUC — unsupervised stage findings

**Status:** unsupervised stage complete; supervised classification not started.
**Dates compared:** 1 Feb – 31 Mar, 2019 and 2024.
**Analysed area:** 7,771.7 km² (97.97% of the dissolved protected-area union).
**Author:** Isaac Mureithi — IMureithi@awf.org

---

## 1. Headline

Between these two dry seasons, **the protected-area network was greener, and
there is no evidence of widespread land-cover loss.** For an audience expecting
a degradation narrative, that is the finding.

**Read "between these two dry seasons" literally.** The annual series (§7.4)
shows the 2019→2024 step is **ordinary** against this landscape's year-to-year
variation — 3 of the 4 comparable five-year steps in the record are at least as
large, and 2020 and 2026 were far greener than either date compared here. This
is a difference between two particular years, not a trend, and never an annual
rate.

The change signal is dominated by a landscape-wide greenness shift between the
two dry seasons, not by conversion and not by fire. Once a change is required to
be large in spectral terms, not explained by fire, and big enough to be a real
patch, the candidate area falls to a small fraction of the landscape.

Three numbers that should **not** be quoted as change, and why:

| Figure | What it actually is |
|---|---|
| 44.6% of area "changed" | k-means cluster membership flipping, mostly at arbitrary cluster boundaries |
| 1,600 km² "substantive change" | pixels moving 2+ cluster ranks — a rank, not a distance |
| 877 km² "residual" | the above, minus fire; inherits the same flaw |

The defensible figure is **263.8 km² (3.39%)** of candidate change — large
spectral movement, no fire to explain it, in patches of 1 ha or more, and
corrected for the observation-depth artefact (§7.5). Candidate, **not
confirmed**: three separate attempts to explain it away have failed (fire,
small missed fires, terrain illumination) and a fourth test could not
adjudicate it, but nothing has yet shown what the greening *is*. Settling that
needs high-resolution imagery.

The uncorrected figure is 277.5 km² (3.57%); quote the corrected one.

**One split that the net figure hides.** About **18.8% of the candidate area —
roughly 52 km² — is browning rather than greening**, and it is concentrated in
particular reserves: Lwala (100% of its candidates), Kidepo Valley (98.8%),
Lopeichubei (84.6%) and Zulia (73.9%), against 0–2% in Timu, Rom and Imatong
(§7.5). Browning is what would indicate clearing or degradation, so those
reserves warrant attention even though the network-level signal is greening.

---

## 2. Processing chain

| Stage | Script | Output |
|---|---|---|
| AOI from WDPA polygons | `R/01_build_aoi.R` | 13 PAs, 7,910.5 km² dissolved |
| HLS coverage screening | `gee/01_hls_coverage_check.js` | clear-observation counts per year and window |
| Burn-season profiling | `gee/01b_burn_season_profile.js` | monthly and per-window burnt fractions |
| MODIS burn rasters | `gee/01c_export_burn_masks.js` | burnt extent and **scar age** on the composite grid |
| Composites | `gee/02_build_composites.js` | 2 × 10-band Int16 GeoTIFF, EPSG:32636, 30 m |
| Composite QC | `R/02_qc_composites.R` | pass/fail verdict, support maps |
| Clustering | `R/03_kmeans_firstlook.R` | k-means classes, change matrix |
| Cluster diagnostics | `R/03b_kmeans_diagnostics.R` | convergence, silhouette, separability |
| k sensitivity | `R/03c_k_sensitivity.R` | k=3 vs k=6, cluster-free estimate |
| Fire attribution | `R/04`, `R/05` | proxy vs MODIS record |
| Change magnitude | `R/06_change_magnitude.R` | change in NDVI rather than cluster ranks |
| Maps | `R/07_maps.R` | eight cartographic figures |
| Sampling-bias test | `R/08_sampling_bias.R` | how much of the greening is observation depth |
| Candidate examination | `R/09_check_candidates.R` | by-PA areas, patch sizes, boundary vs roughness, burn test |
| Depth-adjusted areas | `R/12_depth_adjusted_areas.R` | the headline 263.8 km², and its insensitivity to the coefficient |
| Terrain and observation timing | `gee/01d_export_terrain.js` | elevation, slope, aspect, mean observation day per year |
| Terrain-illumination test | `R/10_terrain_test.R` | whether the candidates are a light artefact |
| Solar-geometry self-check | `R/tests/test_solar_geometry.R` | validates the illumination arithmetic in `R/10` |
| Annual NDVI series | `gee/03_ndvi_timeseries.js` | 11 years × (NDVI, clear looks, burnt) on the composite grid |
| Trend and trajectory analysis | `R/11_ndvi_timeseries.R` | trend vs variability; candidate trajectories |
| Trend-estimator self-check | `R/tests/test_trend_estimators.R` | validates the closed forms and the DiD in `R/11` |

Composites are HLS (Landsat + Sentinel-2 harmonised), Fmask-masked, median
over the window, with NDVI/NDWI/NBR and a clear-observation count carried as
bands. Reflectance and indices are stored ×10000 as Int16 with −9999 nodata.

---

## 3. Why these dates

Both assumptions the project began with turned out to be wrong, and measuring
them changed the choice of window — though not in the direction expected.

**`dry_early` (Dec–Jan) is the clearer window, not the cloudier one.** It
reaches ≥4 clear observations over ~100% of the AOI in 9 of 10 years, against
95–100% in only 6 of 11 for `dry_late`.

**`dry_late` (Feb–Mar) is not peak burning. December is**, then January and
November. Only 1–7% of the AOI burns *during* `dry_late`, against 12–37% during
`dry_early`.

`dry_late` was chosen anyway, for a reason neither assumption anticipated: a
median composite assumes a **stationary surface**. In `dry_early` a pixel's
observations straddle its own fire; in `dry_late` the scars are already down and
ageing uniformly. A stable surface with fewer looks beats a shifting surface
with more.

**2019 and 2024** are the two least-burnt `dry_late` years on record here (19.8%
and 22.9% of the AOI carrying scars), five years apart, both above 96% clear
coverage in each country, with comparable observation density (8.1 and 10.4 mean
clear looks per pixel).

---

## 4. Composite quality

| | 2019 | 2024 |
|---|---|---|
| Mean clear observations | 8.11 | 10.38 |
| Area with <4 observations | 1.92% | 0.40% |
| Median RED / NIR / SWIR1 | 0.095 / 0.198 / 0.240 | 0.087 / 0.191 / 0.230 |
| Median NDVI | 0.324 | 0.333 |

2.03% of the PA union is thinly sampled in one date or the other, and is
excluded from all change statistics. Grids are identical (4559 × 7009, 30 m,
EPSG:32636), verified with `compareGeom` rather than assumed.

**This was tested, not just flagged** — see §7.1. About 28% of the greening is
attributable to 2024's better sampling; the remaining 72% is not.

---

## 5. What the clustering found, and what it cannot support

k-means was fitted once on a pooled sample from both dates and applied to each,
so clusters mean the same thing in both — clustering the dates independently
would have made the change matrix count coincidences.

The clusters are spatially coherent and interpretable: Imatong resolves as
woodland, the Karamoja reserves grade from burnt/bare through dry grassland.

But **the landscape supports about three robust strata, not six**:

| k | Converged | % variance | Mean silhouette | % pixels ambiguous |
|---|---|---|---|---|
| 3 | yes | 71.1 | **0.523** | **16.1** |
| 4 | yes | 78.7 | 0.490 | 18.9 |
| 5 | no | 81.7 | 0.442 | 24.1 |
| 6 | yes | 84.1 | 0.440 | 24.8 |
| 7–10 | no | 86+ | 0.43 falling | 26+ |

Silhouette falls monotonically as k rises while variance explained climbs —
which is why an elbow plot looked reassuring and the structure underneath was
not. At k=6 a quarter of pixels sit ambiguously between clusters; the closest
centroid pair is 1.88 apart in standardised units against 9.04 for the farthest.
Clusters 1–4 span NDVI 0.253–0.390 and are slices of one continuum.

Those boundary pixels flip between dates for no physical reason, which is a
large part of why the raw change figure came out as high as it did.

---

## 6. Fire does not explain the change

Attribution was done against the **MODIS record**, not a proxy. An earlier
NBR-quantile proxy scored kappa 0.086 (2019) and 0.262 (2024) against that
record — at best weak, at worst chance.

| | |
|---|---|
| MODIS burn switched between dates | 2,289 km² |
| Overlap with substantive change | 524 km² |
| Expected under independence | 471 km² |
| **Enrichment** | **1.11×** |
| Risk ratio | 1.17× |

Burning barely predicts where the clusters changed. Burnt ground covers ~20% of
the AOI in both dates, but it is largely the *same* ground in both, and fire is
therefore a constant rather than a driver for those pixels.

### Why burn scars are nearly invisible here

| Days since burn | Median NBR |
|---|---|
| 0–20 | **+0.294** |
| 40–60 | +0.088 |
| 80–100 | −0.020 |
| 120–160 | −0.066 |
| *(unburnt ground)* | *−0.001* |

The **youngest** scars look greenest. This is not scar recovery — it is
compositing. A pixel that burnt in mid-March has most of its February–March
observations *predating its own fire*, so the median shows unburnt ground. A
pixel that burnt in November has every observation after the fire. A Feb–Mar
median therefore hides recent fires and shows only old, faded ones.

This is a methodological finding worth carrying into any future dry-season
composite work in fire-driven systems: **a seasonal median composite is a poor
instrument for detecting fire that occurred within or just before its window.**

---

## 7. The change, measured properly

Cluster "steps" are a rank, not a distance. The k=6 centroids sit at NDVI 0.252,
0.259, 0.307, 0.387, 0.521, 0.792 — so a 1→3 move counts as two steps while
spanning 0.055 NDVI, and 5→6 counts as one while spanning 0.271.

| Cluster steps moved | % of area | Median \|ΔNDVI\| |
|---|---|---|
| 0 | 55.6 | 0.033 |
| 1 | 23.9 | 0.045 |
| 2 | 14.2 | 0.060 |
| 3 | 6.0 | 0.072 |
| 4+ | 0.37 | **0.220** |

Only the 0.37% moving four or more steps shift substantially.

### Sensitivity to k

| | k=3 | k=6 | ratio |
|---|---|---|---|
| Changed cluster | 1,568 km² (20.2%) | 3,467 km² (44.6%) | 2.21× |
| …not explained by fire | 1,076 km² | 2,289 km² | 2.13× |
| …and \|ΔNDVI\| > 0.15, ≥1 ha | 146.8 km² | 226.2 km² | **1.54×** |

Filtering removes most of the k-dependence but not all of it. A 50% swing from a
parameter the data cannot justify is not a measurement, which is why the
headline estimate below uses **no clustering at all**.

### The cluster-free estimate

Criteria, all defined in sensor units and reproducible without tuning:
|ΔNDVI| > 0.15, no MODIS burn-state switch to explain it, and contiguous
patches ≥ 1 ha.

> ### 277.5 km² — **3.57%** of the analysed area
>
> in contiguous patches of 1 ha or more.
>
> **Corrected for observation depth in §7.5: 263.8 km² (3.39%).** That is the
> figure to quote; the one above is the uncorrected input to it.

This replaces both cluster-dependent versions (146.8 km² at k=3, 226.2 km² at
k=6). It is **larger** than either, because it drops the requirement that a
pixel also change cluster — the criterion is now purely spectral, so it catches
large spectral moves that happened to stay inside one cluster. That is the
intended behaviour: the clustering was never the thing doing the measuring.

A note against over-reading it: the layer is size-filtered by construction, so
its spatial coherence ratio (20.9, with 1.3% isolated pixels) is **not**
independent evidence of anything. A layer built from ≥1 ha patches is coherent
because it was made that way.

### Direction of change

| | Median ΔNDVI | % greening |
|---|---|---|
| Background | 0.012 | 59.2 |
| Residual (fire removed) | 0.037 | 65.8 |

Landscape-wide median ΔNDVI is **+0.0137** — 2024 greener than 2019.

Land-cover conversion produces clearing as well as regrowth, so a conversion
signal would sit near 50% greening. A landscape moving one way together is a
wetter or later-senescing season, and no burn mask can remove it.

**But do not read that as "nothing browned".** Among the *candidate* patches
specifically, **18.8% are browning**, and §7.5 shows they are not spread evenly:
four reserves carry candidates that are predominantly or entirely browning while
the three that drive the greening are almost purely greening. The
landscape-level direction and the per-reserve direction are different findings.

### 7.1 Is the greening just unequal sampling?

2019 carries 8.2 mean clear looks against 2024's 10.5. If cloud screening
removed proportionally more of the greener, cloudier days from 2019, that alone
would manufacture the result. `R/08_sampling_bias.R` tests it three ways.

The decisive test regresses ΔNDVI on the per-pixel **depth gap** with
protected-area fixed effects. Differencing the same pixel across two dates
cancels everything time-invariant about it — terrain, soil, land cover — leaving
the depth gap as what varies:

| | ΔNDVI |
|---|---|
| Mean observed | +0.0176 |
| Attributable to the depth gap | **+0.0049 (28%)** |
| Remaining at zero depth gap | **+0.0127 (72%)** |

**The greening is real but was overstated.** Quote roughly **+0.013**, not
+0.018, and to one significant figure — estimates range from +0.009 to +0.019
across balanced subsets.

Two methodological notes, because the tests disagreed and the disagreement is
instructive:

- A *within-date* regression of NDVI on observation count comes out **negative**
  in both years, which would suggest the mechanism cannot operate at all. That
  comparison is between different pixels, and the ones with fewer looks are the
  cloudy highlands, which are also the green ones. It measures geography as much
  as sampling, and is the weaker estimate.
- Restricting to "balanced" subsets does **not** hold the landscape fixed: mean
  looks climb from 8.2/10.5 to 14.3/14.5 as the rule tightens, so the strictest
  subset is the 7.4% of the landscape with the clearest skies. Those subsets
  confound depth balancing with spatial selection, which is why their medians
  bounce rather than trend. They establish only that every subset still greens.

---

### 7.2 The candidate patches, examined

`R/09_check_candidates.R`. Four findings, one of which overturned my own
prediction.

**They concentrate in the small hill reserves, not the flagship parks:**

| Protected area | Area km² | Candidate km² | % of PA |
|---|---|---|---|
| Timu | 122.1 | 32.3 | **26.5** |
| Rom | 108.9 | 25.6 | **23.5** |
| Lwala | 58.7 | 6.1 | 10.3 |
| Imatong | 1001.2 | 72.1 | 7.2 |
| Kidepo (SS) | 2828.0 | 59.2 | 2.1 |
| Kidepo Valley | 889.6 | 4.2 | **0.5** |

**They sit nearer protected-area boundaries** — median 1,700 m against 3,026 m
for background. In a logistic model with standardised predictors, distance to
boundary (−0.85) outweighs local spectral roughness (+0.49). That matters
because a legal boundary has no inherent spectral gradient, whereas
misregistration produces artefacts at *image* edges. So boundary proximity is
the stronger pattern — though roughness is far from negligible, and boundaries
often follow rivers and ridges, which are image edges too.

**They green, strongly and almost exclusively** — median ΔNDVI **+0.18** with
**82%** greening, against +0.012 and 59% for background. Clearing would brown;
conversion would do both.

**They are not recovering fires.** This was my leading hypothesis — MODIS at
500 m cannot see a few-hectare fire, so small 2019 burns recovering by 2024
would look exactly like this. It is wrong: candidates are only **0.73×** as
likely as background to have been in the 2019 burnt tail, and their 2019 NBR
(+0.081) was already *above* the landscape median (−0.009). They were greener
than average before and greener still after.

**What remains, none of it established:** genuine woody thickening or regrowth;
differential terrain illumination, since candidates favour rough, hilly ground
in reserves that are mostly forested massifs; or real land-use change. The
illumination explanation is the one that can be settled from orbit, and §7.3
sets out the test. Separating the other two needs high-resolution imagery or a
multi-year series. One caution on the patch statistics: eight-connected
labelling threads scattered clumps into one component, so the "largest patch,
3,064 ha" is a sprawl across a hillside, not a block. The median patch is
2.34 ha.

---

### 7.3 Terrain illumination is not the explanation — settled

**Status: run.** `R/10_terrain_test.R` against the NASADEM terrain export. The
expectations below were recorded *before* the export existed, and the result
matched them.

**The verdict: illumination explains +1.3% of the greening, and 0% of the
candidate area is illumination-explicable. The candidate set survives intact.**
(This test ran against the uncorrected 277.5 km² raster; the depth correction in
§7.5 is independent of it and reduces the area to 263.8 km² for unrelated
reasons.)

**The decisive argument needs no regression.** Turn the question round: how
strongly would NDVI have to respond to illumination for illumination to explain
the observed +0.0177 mean ΔNDVI? That is just the greening divided by the
measured shift in cos(i):

| | Mean Δcos(i) | Required coefficient |
|---|---|---|
| Computed angles | +0.0024 | **+7.50** NDVI per unit cos(i) |
| Recorded angles | −0.0065 | **−2.73** |

The measured response is **+0.099** (differenced) or **−0.119**
(cross-sectional), and NDVI spans about one unit end to end — so the required
value is 63× and 23× the largest response this data shows, and beyond what is
physically available at all. This rests on arithmetic, not on any regression's
identifying assumptions.

**The mechanism was never available.** Both composites use the same calendar
window, so the sun sits in nearly the same place in each. Measured: observation
timing shifted **+2.0 days** between 2019 and 2024, moving the mean solar zenith
by **−0.29°** (computed) or **+0.79°** (recorded) — both under a degree, and the
two sources need not agree on the sign of something that small. cos(i) moved by
0.24% of its range.

**No aspect fingerprint.** Candidate enrichment across the eight compass octants
spans just **1.15×** (strongest W at 1.04×), and the sun-facing and shaded groups
both come out at 1.00×. Illumination has a specific signature — a preference for
lit or shaded slopes — and it is absent.

**What the candidates *are* is higher and steeper**, confirming the hill-reserve
pattern from §7.2 without supporting an artefact:

| | Candidates | Background |
|---|---|---|
| Median slope | 8.0° | 4.0° |
| Mean slope | 12.6° | 8.1° |
| Share above 10° | 45.0% | 25.6% |
| Median elevation | 1,406 m | 1,151 m |

In a *differenced* regression, slope predicting change is a statement about
where the vegetation is, not about the light.

**Route 2 is refuted too, and in the opposite direction to the hypothesis.** The
idea was that Fmask cannot tell terrain shadow from cloud shadow, so the §7.1
depth artefact should concentrate on steep ground. It does the reverse: the
depth effect is **+0.0029 per look on flat ground against +0.0009 on steep**, a
ratio of **0.32×** (interaction p < 2e-16). The depth confound is real but it
lives on flat ground, so it needs no separate terrain correction.

A useful by-product: setting both the depth gap and slope to zero puts the
greening at **+0.0113**, against §7.1's +0.0127 from a different specification.
Two routes landing within 0.0014 of each other is a genuine cross-check.

**The hypothesis as originally stated:** reflectance on a slope depends on how
the sun strikes it; the candidates favour rough, hilly ground, so some of the
"change" could be light rather than land.

**Why the hypothesis is already in trouble, on mechanism.** Both composites use
the *same calendar window*, 1 Feb – 31 Mar. The sun is therefore in nearly the
same place in both, so the local illumination geometry cos(i) is nearly the same
in both, so it very largely **cancels in the difference**. Illumination sets the
*level* of reflectance on a slope; it cannot by itself produce a *change*
between two same-season composites. For it to do so, something about the
observations must differ between the years, which leaves two routes:

- **Route 1 — timing.** If one year's clear looks fall earlier or later in the
  window, the sun sat at a different declination and cos(i) shifts, most on
  steep slopes. The export's `obs_doy_<year>` bands measure this per pixel.
- **Route 2 — shadow masking.** Fmask cannot tell terrain shadow from cloud
  shadow, so shaded slopes lose observations. §7.1 already attributed 28% of the
  greening to observation depth; if that artefact concentrates on steep ground,
  terrain is its mechanism. This is the version most likely to survive — and it
  is not really about illumination at all, but about where the §7.1 confound
  lives.

**What the test can see.** Validated by simulation before the export, in
`R/tests/test_solar_geometry.R`: the differenced estimator recovered a planted
coefficient of 0.30 as **0.3002** over 40 replications with no detectable bias
(t = +0.09), and its reported standard error matched the empirical spread. At
n = 400,000 the standard error on the coefficient is about 0.007, so the test
has ample power.

**What it cannot see is a large effect, because there isn't one available.** The
greening illumination can explain is the coefficient times the mean shift in
cos(i), and for a timing shift of a few days that shift is of order 0.004:

| Coefficient | Greening explained | Share of +0.0127 |
|---|---|---|
| 0.1 | +0.0004 | 3% |
| 0.3 | +0.0012 | 9% |
| 0.5 | +0.0020 | 16% |
| 1.0 | +0.0040 | 32% |

NDVI's cross-sectional spread across aspects in rugged terrain is rarely more
than ~0.1 per unit cos(i), so coefficients near 1.0 are not physically
available. **The ceiling is single-digit percent unless the two years were
observed many days apart** — which is exactly what Phase 2 of the script
measures, and why the obs_doy bands matter more here than the DEM does.

**The diagnostic to watch is aspect, not slope.** "Candidates are steeper" is
ambiguous: woody vegetation genuinely prefers hills in Karamoja. But there is no
ecological reason for regrowth to prefer west-facing over east-facing ground at
the same slope and elevation, whereas illumination has exactly that signature —
at a mid-morning overpass the sun sits east-southeast (computed: zenith 26–34°,
azimuth 89–130° across the window), so east-facing slopes are lit and
west-facing ones are not. A north–south asymmetry would *not* be illumination:
at 3–4°N the sun passes nearly overhead and both are lit about equally.

**What this means for the supervised stage.** A topographic correction
(C-correction, SCS+C) is **not** needed before the Random Forest. That was the
contingency if illumination had proved material; it did not.

**What remains open.** Two explanations this pipeline cannot separate: genuine
woody thickening or regrowth, and real land-use change. Both would need
high-resolution imagery for a sample of patches. §7.4 asks whether the annual
series can help.

**One caution on the recorded angles.** The first run reported recorded and
computed zenith differing by 20.1°, which would have discredited the
computation. The computation was right; the reference was contaminated by
antimeridian scenes leaking through `filterBounds()` — see §10. Restricted to
the correct MGRS zone, the two agree within 0.6°.

---

### 7.4 The annual series: the step is ordinary, and the candidates are not singled out

**Status: run.** `R/11_ndvi_timeseries.R` over the 2016–2026 annual dry-season
series (NDVI, clear looks and burnt mask, three 11-band stacks).

**Two results, one about the landscape and one about the candidates.**

#### The 2019→2024 step is ordinary

On a **fixed** set of 38,622 pixels — the same pixels every year, so composition
cannot move the series — the annual means are:

| 2016 | 2017 | 2018 | 2019 | 2020 | 2021 | 2022 | 2023 | 2024 | 2025 | 2026 |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.365 | 0.322 | 0.366 | **0.366** | **0.481** | 0.406 | 0.349 | 0.348 | **0.387** | 0.372 | **0.511** |

- The headline step, 2019→2024, is **+0.0210**.
- The year-to-year standard deviation is **0.0761**; the largest single-year
  swing is **+0.139**.
- **3 of the 4 comparable five-year steps in this record are at least as
  large.**
- The series spans 0.322 to 0.511 — a range about **nine times** the step.

2020 and 2026 are genuinely far greener than the two chosen years, and that is
not a coverage artefact: the composition effect (fixed minus varying footprint)
is at most **0.009** inside the panel. Two dates five years apart in this
landscape cannot establish a trend, and the figure must never be quoted as an
annual rate.

#### The candidates are not distinguishable over the series

This is where the analysis nearly went wrong, and the correction is the most
important methodological point in this document.

The candidates were **defined** as pixels with |ΔNDVI| > 0.15 between 2019 and
2024. Selecting pixels on a difference between two years guarantees that they
will look low in the first year, high in the second, and rising in between —
whether or not anything real is happening. Fitting a trend across a window that
*contains* those years therefore partly restates the selection.

Measured in-sample, the candidates looked convincing: a trend **+0.00517 ±
0.00023 NDVI/yr** faster than their background with PA fixed effects — 22
standard errors — which the script's own text called "the clean number". The
tell that something was wrong was the companion figure: an excess 2019 dip of
−0.067, **three times the headline step it was supposed to explain.**

So the same rule was applied to other year pairs five years apart, as a placebo.
Those sets were not selected on 2019 or 2024, so whatever they show is what the
rule manufactures from a landscape where nothing has been singled out:

| Pair | Pixels flagged | % flagged | Leave-out trend gap | Excess dip at first year |
|---|---|---|---|---|
| 2018→2023 | 258 | 0.67 | **+0.00582** | −0.0845 |
| **2019→2024 (real)** | 2,041 | 5.28 | **+0.00424** | **−0.0845** |
| 2020→2025 | 64 | 0.17 | **+0.01376** | −0.1217 |
| 2021→2026 | 11,310 | 29.28 | +0.00034 | −0.0271 |

**The real pair sits inside the placebo range on both measures, and below the
placebo mean on the trend (+0.00424 against +0.00664).** Its excess dip is
*identical* to the 2018→2023 placebo. Both apparent findings belong to the
selection rule, not to the ground.

**Conclusion: the series adds no independent support for the candidate patches.**
They still rest on the two-date comparison alone. It does not refute them
either — it says this instrument cannot see them.

One caveat on the placebo, which strengthens rather than weakens it. The flagged
fractions differ enormously because the year pairs differ in their
landscape-wide NDVI change: 2021→2026 rose 0.105 across the whole landscape, so
a fixed 0.15 threshold catches 29% of it and is barely selective at all. The
comparisons worth weighing are therefore with the pairs of *similar or greater*
selectivity — and both of those, 2018→2023 at 0.67% and 2020→2025 at 0.17%,
produce a **larger** trend gap than the real pair. Only the barely-selective
pair gives a smaller one.

The probe uses a simplified one-sided rule (ΔNDVI > 0.15, no fire switch, no
minimum mapping unit) so every pair is treated identically; that is why it flags
5.28% for the real pair where §7's full rule gives 3.57%.

#### Why no absolute trend is quoted, and a correction to an earlier claim

Observation depth rises across the record (correlation with year **+0.64**:
Sentinel-2B from 2017, Landsat 9 from 2021), and §7.1 established that *within a
pixel* more looks produce higher NDVI (+0.0033 per look). Before the data
existed, this document predicted that a denser archive alone could manufacture
the whole signal.

The data says something more awkward. Across years, depth correlates with mean
NDVI at **−0.52** — the wrong sign. Wet years are cloudy years: they yield
*fewer* clear looks and *greener* vegetation, and that mechanism dominates the
annual means. So there are two confounds of **opposite sign** — a within-pixel
compositing artefact pushing a trend up, and a rainfall–cloud relationship
pushing it down — and neither can be separated from the other or from a real
trend with nine annual points. The earlier claim that depth would inflate the
trend is withdrawn; its net direction is unknown.

The same problem defeats correcting the trend map. Estimated within year across
pixels, the depth coefficient comes out **−0.00049 per look**, opposite in sign
to the differenced +0.0033 — exactly the failure §7.1 warned of, where a
cross-pixel slope measures which ground gets seen rather than the artefact.
Applying it would have *inflated* the trend (+0.0059 to +0.0064 per year). **No
depth correction is applied**, and the map is published raw.

#### The trend map, and what it covers

`outputs/timeseries/ndvi_trend_slope.tif`, mean slope **+0.0059 NDVI/yr** over
the balanced panel — a level that is not interpretable, for the reasons above;
read the pattern only. Requiring all nine years to be adequately supported keeps
just **19.5% of the AOI**, and that fifth is flatter and less cloudy than the
whole landscape.

#### Fire

Cumulative burnt fraction runs 19.7% (2019) to 43.9% (2018), but its correlation
with year is only **+0.04** — fire varies between years without trending, so it
does not bias the series direction.

---

### 7.5 The candidate areas, depth-adjusted

`R/12_depth_adjusted_areas.R` closes the oldest outstanding caveat in this
document: §7.1 corrected the headline ΔNDVI for observation depth but never
corrected the *areas*, which were left flagged as "likely a little generous".

**The rule is reproduced before it is adjusted.** The script rebuilds §7's
cluster-free candidate set from scratch and gets **277.52 km²** against the
published 277.5 — a difference of 0.02 km². An adjusted area means nothing
unless the unadjusted one can be recovered first, so the script stops if it
cannot.

**The adjustment is per-pixel and slope-dependent**, using §7.3's finding that
the depth response is about three times stronger on flat ground:

`ΔNDVI.adj = ΔNDVI − max(b_gap + b_inter · slope, 0) · obs.gap`

The slope *main* effect is deliberately not subtracted — in a differenced
regression that is a statement about where the vegetation is, not an artefact.
The coefficient is clamped at zero above its 29.3° break-even, affecting 4.7% of
pixels, because a negative depth response has no mechanism.

**Result: 277.5 km² → 263.8 km², a fall of 13.7 km² (−4.9%).** Of that, 25.5 km²
was removed and 11.8 km² added.

**And it is not sensitive to the coefficient**, which was predicted in advance:
the adjustment is of order 0.006 NDVI against a threshold of 0.15, thirty times
smaller than the bar it has to move pixels across.

| Adjustment strength | Area | Change |
|---|---|---|
| ×0 (check) | 277.52 km² | 0.0% |
| **×1 (estimated)** | **263.80 km²** | **−4.9%** |
| ×2 | 254.69 km² | −8.2% |
| ×3 | 251.54 km² | −9.4% |

A threefold overcorrection moves the area by under 10%. The caveat was real but
small.

**A cross-check worth recording.** Re-estimating §7.1's simpler specification on
a fresh sample gives **+0.002120 NDVI per clear look** against its own implied
~0.00214 — essentially identical. With the slope interaction the attributable
share is 32% and the corrected greening **+0.0120**, against §7.1's 28% and
+0.0127. Two specifications, two samples, agreeing within 0.0007.

#### Why some reserves gain area, and the split this uncovered

The adjustment does not shrink every reserve. Its direction obeys

`sign(area change) = −sign(observation gap) × sign(ΔNDVI)`

because the threshold is **two-sided**. A reserve seen *less* often in 2024 has
its 2024 NDVI biased downward, so a **greening** pixel there is pushed further
past +0.15 and gained — while a **browning** pixel in the same reserve is pushed
back towards zero and lost.

Both terms are needed. Taking the gap alone, three of thirteen reserves look
like they break the rule; with the browning share beside it, all thirteen follow
it:

| | Clear looks 2019 → 2024 | Browning share of candidates | Area change |
|---|---|---|---|
| Timu | 20.6 → 14.8 (**−5.8**) | 0.0% | **+18.4%** |
| Kidepo Valley | (+1.6 among candidates) | **98.8%** | **+24.6%** |
| Lwala | 9.6 → 6.3 (−3.5) | **100.0%** | −6.9% |
| Lopeichubei | (−0.8) | **84.6%** | −0.9% |
| Imatong | 6.2 → 10.3 (+5.5) | 2.2% | −12.8% |
| Kidepo | 6.8 → 10.5 (+3.7) | 25.1% | −10.2% |

The three rows that the gap alone could not explain carry **94.5% browning
candidates against 12.2% elsewhere.** Two earlier explanations for them — the
slope clamp, then a mismatch between whole-reserve and candidate-pixel gaps —
were both tested and both wrong; this one accounts for every row.

**And it exposes something that matters more than the adjustment did.** About
**18.8% of the candidate area is *browning*, not greening — roughly 52 km²** —
and it is heavily concentrated:

| Reserve | Candidate area | Browning share | Browning area |
|---|---|---|---|
| Zulia | 21.2 km² | 73.9% | ~15.7 km² |
| Kidepo | 59.2 km² | 25.1% | ~14.9 km² |
| Lopeichubei | 10.1 km² | 84.6% | ~8.5 km² |
| Lwala | 6.1 km² | 100.0% | ~6.1 km² |
| Kidepo Valley | 4.2 km² | 98.8% | ~4.1 km² |
| Timu, Rom, Imatong, Nyangea, Agoro | 146.5 km² | 0–2.2% | ~2.3 km² |

So the landscape-level finding — greener overall, no evidence of widespread loss
— is true and **conceals a split**. The large hill reserves that drive the
greening (Timu, Rom, Imatong) are almost purely greening, while Lwala, Kidepo
Valley and Lopeichubei have candidate change that is almost entirely the other
way. For management the split is the more useful fact, and **browning candidates
are the ones that would indicate clearing or degradation**, so they deserve the
first look when high-resolution imagery is obtained (§11).

Two caveats carry over: these are still *candidates*, not confirmed change; and
the reserves concerned are small, so a few patches crossing the
minimum-mapping-unit threshold move their percentages a long way.

---

## 8. Maps and figures

All generated by `R/07_maps.R` into `outputs/maps/` (gitignored — re-run the
script to regenerate). Built with ggplot2 + tidyterra, ggspatial for scale bars
and north arrows, patchwork for panels, rnaturalearth for country context.

| File | Shows |
|---|---|
| `map_01_study_area.png` | The 13 protected areas by country, with context and labels |
| `map_02_composites.png` | False-colour SWIR2/NIR/GREEN composites, 2019 vs 2024 |
| `map_03_observation_support.png` | Clear looks per pixel — the honesty map |
| `map_04_kmeans_classes.png` | k-means classes both dates, semantic palette |
| `map_05_fire_state.png` | MODIS burnt area by date |
| `map_06_scar_age.png` | Days since burn at the window's end |
| `map_07_dndvi.png` | Greenness difference, diverging scale |
| `map_08_candidate_change.png` | Candidate change — the **uncorrected** 277.5 km² set, since `07_maps.R` predates §7.5. Re-run it against `change_strong_clusterfree_depthadj.tif` before using this map in a report |

Zoom panels of the largest candidate patches, before and after, are written by
`R/09_check_candidates.R` into `outputs/candidates/`.

**Terrain test figures** — `R/10_terrain_test.R` into `outputs/terrain/`:

| File | Shows |
|---|---|
| `terrain_aspect_rose.png` | Candidate enrichment by compass octant — the diagnostic figure for §7.3. A flat rose is the absence of an illumination fingerprint |
| `terrain_candidates_over_hillshade.png` | Candidates over a hillshade of the worst-affected reserve, lit from the east-southeast as at the overpass. Candidates clinging to one side of every ridge would have meant illumination |
| `terrain_slope_elevation_plane.png` | Candidate frequency across slope and elevation — what the candidates genuinely are |
| `terrain_dndvi_vs_dcosi.png` | ΔNDVI against the change in illumination; the null made visible |

**Annual series figures** — `R/11_ndvi_timeseries.R` into `outputs/timeseries/`:

| File | Shows |
|---|---|
| `ts_01_annual_series.png` | Mean NDVI by year on a fixed pixel set, with the two headline years circled, point size for observation depth and colour for fire. The red dashed line is the varying-footprint series |
| `ts_02_candidate_trajectories.png` | Candidate against background year by year |
| `ts_03_slope_distributions.png` | Per-pixel trend distributions for the two groups |
| `ts_04_trend_map.png` | Uncorrected NDVI trend 2018–2026 over the balanced panel (19.5% of the AOI) |

**Depth adjustment** — `R/12_depth_adjusted_areas.R` writes
`outputs/qc/depthadj_candidate_change.png`, showing which candidate pixels the
observation-depth correction removes and which it adds.

Read `ts_02` and `ts_03` with §7.4's selection-bias warning attached: they show
a real difference between the two groups, but the placebo test establishes that
a difference of that size is what the selection rule produces on its own.

Diagnostic plots from earlier scripts remain in `outputs/qc/` and
`outputs/classification/`. Those are working figures, not presentation ones —
no scale bars, no north arrows — and should not go in a report.

---

## 9. What this does and does not establish

**Established:**

- The two composites are sound, matched and comparable, verified rather than assumed.
- Fire does not drive the change signal, tested against recorded fires.
- The raw unsupervised change figures are inflated by cluster-boundary effects.
- The dominant signal is a modest, one-directional greening.
- Terrain illumination does not explain it (§7.3): the response required is
  23–63× what is measured, and beyond what is physically available.
- The 2019→2024 step is **ordinary** against this landscape's year-to-year
  variation (§7.4), so it is a difference between two particular years and not
  a rate.

**Not established:**

- That the remaining candidate patches *are* land-cover change. They are
  candidates: a large greenness change with no *recorded* fire to explain it.
  §7.2 rules out recovering small fires and shows the patches green from an
  already-vegetated baseline, concentrated in the small hill reserves; §7.3
  rules out terrain illumination; §7.4 finds the annual series cannot
  distinguish them from what the selection rule produces on its own. What that
  greening *is* remains open, and three tests have now failed to settle it.
- Any thematic class labels. k-means finds spectral clusters; the tentative
  names in the scripts are read off centroid signatures and are unverified.
- Anything about causes. Greening is consistent with rainfall differences, woody
  encroachment, grazing changes or recovery, and this analysis separates none
  of them.

**Carry these caveats:**

1. 2019 has ~25% fewer clear observations than 2024. Tested (§7.1): it inflates
   the greening by 28–32%, so the corrected figure is **+0.012 to +0.013 ΔNDVI**
   depending on specification. The candidate areas have now been recomputed on a
   depth-adjusted basis (§7.5): **263.8 km² against 277.5 km², −4.9%.** This
   caveat is closed.
2. MODIS burned area is 500 m information — good for landscape attribution, not
   for per-patch verdicts.
3. **Two dates cannot separate inter-annual variability from trend, and the
   series confirms it rather than fixing it** (§7.4). The step is within one
   year-to-year standard deviation and 3 of 4 comparable steps are at least as
   large. No absolute trend is recoverable: two confounds of opposite sign
   (within-pixel compositing, and wet-years-are-cloudy) cannot be separated with
   nine annual points.
4. Protected-area polygons overlap by 584 km², so per-PA areas double-count
   unless overlaps are assigned. The §7.2 per-PA table inherits this.
5. **Anything measured on the candidate set is subject to selection bias**
   (§7.4). They were defined by a 2019→2024 difference, so they will look low in
   2019, high in 2024 and rising in between for free. Any future claim about
   them needs a placebo or a leave-out design; an in-sample contrast on these
   pixels is not evidence, however many standard errors it carries.

---

## 10. Data quality problems found and fixed

Recorded because each was invisible until specifically tested for.

| Problem | How it surfaced | Fix |
|---|---|---|
| All six reflectance bands exported as zero | QC printed values; a range check had passed because zero is "plausible" reflectance | HLS v002 in GEE is already reflectance-scaled; the 0.0001 factor was applied twice. Band-mean smoke test added before export, and a degeneracy check (`min == max`) in QC |
| Band means 5× too low | Median over valid pixels disagreed with the mean | GEE writes pixels outside the clip as 0, not the sentinel; ~70% of the raster was zero-fill. All statistics now computed AOI-masked; GEE clips before unmasking |
| MODIS binary burn bands exported as zero | Age band was correct from the same collection and dates | Likely an `ee.Algorithms.If` branch-projection mismatch. Binary bands now derived from the age band |
| k never actually chosen | The script's own comment said to revisit it after the elbow; it was never acted on | `03b` supplies convergence, silhouette and separability; k=6 retained as reconnaissance only |
| **Antimeridian scenes leaking into every scene query** | `R/10` reported recorded and computed solar zenith differing by 20.1° — a gap big enough to discredit the computation | 235 of 377 scenes in the sun-geometry export were MGRS **zone 01 and zone 60** tiles, acquired at a mean of 18:00 UTC against 07:24 for the real ones. A tile straddling 180° has a footprint that wraps in longitude, so `filterBounds()` matches it against any geometry anywhere. All 235 were Sentinel-2; every Landsat scene was genuine. **No pixel-level result is affected** — a leaked scene has no data over Kidepo, so it contributes a masked pixel and drops out of every median, count and mean. What it corrupts is anything computed from per-*scene* metadata. Restricting to `T36` tiles brought recorded and computed zenith to within 0.6° |
| **`n_scenes` columns are inflated** | Same cause as above | The scene counts in `kidepo_hls_coverage_*.csv` (§3) and `kidepo_ndvi_annual_aoi.csv` include antimeridian tiles — e.g. dry_late 2019 reports 141 Sentinel-2 scenes where only ~43 are over the AOI. The window and year choice rested on `mean_obs` and `pct_ok`, which are per-pixel and unaffected, so no decision changes. Do not quote the scene counts |
| **Zero-fill recurred, and the earlier fix was mis-described** | `R/10`'s own range assertion refused to run: elevation minimum 0 in a landscape whose valley floor is at 709 m | The 2024 note above said clipping before unmasking fixed this. It does not. Earth Engine exports the bounding **rectangle** of a non-rectangular region, fills outside the region with 0, and writes **no nodata tag** — 70.5% of the terrain export. That happens in the export, downstream of the image, so no ordering of `clip` and `unmask` prevents it. The real and only fix is the one already in use: AOI-mask before computing. `R/10` and `R/11` now do that *before* their encoding checks, and the misleading comments in the GEE scripts are corrected |
| **Annual means compared across different sets of pixels** | The exported annual table read 0.491 NDVI for 2020 and 0.528 for 2026 against ~0.38 for every other year — and both are low-coverage years | A year's mean over whatever pixels happened to be clear is a mean over a *different landscape* each year; in 2020 only 75.9% of the AOI was adequately covered, and the surviving part is not a random part. `R/11` now computes the series over a **fixed** pixel set and prints the varying-footprint version beside it as a `composition` column |

---

## 11. Next

The unsupervised stage is finished. Everything that could be tested from orbit
has been: fire (§6), sampling depth (§7.1), the patches themselves (§7.2),
terrain illumination (§7.3) and the annual series (§7.4). Three of those were
attempts to explain the candidate change away, and none succeeded; the fourth
showed the series cannot adjudicate it at all.

1. **Check candidate patches against high-resolution imagery.** This is now the
   only remaining way to say what the change *is*. Nothing in this pipeline can
   do it — it needs eyes on imagery the project does not carry. A sample of
   20–30 patches stratified by protected area would settle whether they are
   woody thickening, regrowth or land-use change.
   **Start with the browning patches** (§7.5): ~52 km², concentrated in Lwala,
   Kidepo Valley, Lopeichubei and Zulia. They are a fifth of the candidate area
   but they are the fifth that could indicate clearing or degradation, so they
   carry the most management consequence per patch inspected.
2. **Agree a class scheme** before labelling. The k=3 strata are a sensible
   skeleton; burnt vs bare and grassland vs wooded grassland need training data.
3. **Training data** for the supervised stage — the open decision. Note what
   §7.3 settled: no topographic correction is needed first.
4. **If the candidate set is used in the supervised stage at all**, carry §9
   caveat 5 with it. Those pixels were chosen on a two-year difference, so any
   statistic computed on them in-sample is contaminated by that choice.
